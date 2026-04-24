from ast import literal_eval
import logging
import re
from typing import Dict

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph
from multiqc.utils import mqc_colour

log = logging.getLogger(__name__)

PATHOGEN_MAP = {
    "2697049": "SARS-CoV-2",
    "208893": "RSVa",
    "208895": "RSVb",
}


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Freyja",
            anchor="freyja",
            href="https://github.com/andersen-lab/Freyja",
            info="Recovers relative lineage abundances from mixed samples.",
            extra="""
            Freyja is a tool to recover relative lineage abundances from mixed samples from a
            sequencing dataset and uses lineage-determining mutational "barcodes" derived from the UShER global
            phylogenetic tree to solve the constrained (unit sum, non-negative) de-mixing problem.
            """,
            doi="10.1038/s41586-022-05049-6",
        )

        # data_by_taxid[taxid][s_name] = {lineage: abundance}
        data_by_taxid: Dict[str, Dict[str, Dict[str, float]]] = {}

        for f in self.find_log_files("freyja", filehandles=True):
            s_name = f["s_name"]
            taxid_match = re.search(r"-(\d+)$", s_name)
            if not taxid_match:
                raise ValueError(f"Could not extract taxid from s_name '{s_name}'")
            taxid = taxid_match.group(1)
            if taxid not in PATHOGEN_MAP:
                raise ValueError(f"Unknown taxid '{taxid}' in s_name '{s_name}'")

            # Strip the "-<taxid>" suffix from the sample name
            s_name = re.sub(r"-\d+$", "", s_name)

            # This will not raise because the search pattern requires it to be there:
            summarized_line = next(line for line in f["f"] if line.startswith("summarized\t"))
            dict_str = summarized_line.split("\t")[1].strip()
            try:
                sample_dict: Dict[str, float] = dict(literal_eval(dict_str))
            except (ValueError, SyntaxError):
                log.error(f"Error parsing 'summarized' line for '{s_name}': {dict_str}, skipping sample")
                continue
            if not sample_dict:
                log.debug(f"No data in the 'summarized' line for '{s_name}': {dict_str}, skipping sample")
                continue

            if taxid not in data_by_taxid:
                data_by_taxid[taxid] = {}
            if s_name in data_by_taxid[taxid]:
                log.debug(f"Duplicate sample name found for taxid {taxid}! Overwriting: {s_name}")
            data_by_taxid[taxid][s_name] = sample_dict
            self.add_data_source(f, s_name)

        # Remove filtered samples
        for taxid in list(data_by_taxid.keys()):
            data_by_taxid[taxid] = self.ignore_samples(data_by_taxid[taxid])
            if not data_by_taxid[taxid]:
                del data_by_taxid[taxid]

        if not data_by_taxid:
            raise ModuleNoSamplesFound

        total_samples = sum(len(v) for v in data_by_taxid.values())
        log.info(f"Found {total_samples} reports")

        self.add_software_version(None)

        # Sort each organism's samples for reproducible order
        for taxid in data_by_taxid:
            data_by_taxid[taxid] = dict(sorted(data_by_taxid[taxid].items()))

        self.scale = mqc_colour.mqc_colour_scale("plot_defaults")

        for taxid in sorted(data_by_taxid.keys()):
            organism = PATHOGEN_MAP[taxid]
            samples = data_by_taxid[taxid]

            top_lineages_dict: Dict[str, Dict] = {}
            all_lineages: list = []
            seen_lineages: set = set()

            for s_name, sample_data in samples.items():
                top_lineage, top_lineage_value = max(sample_data.items(), key=lambda xv: xv[1])
                top_lineages_dict[s_name] = {
                    f"Top_lineage_freyja_{taxid}": top_lineage,
                    f"Top_lineage_freyja_{taxid}_percentage": top_lineage_value,
                }
                if top_lineage not in seen_lineages:
                    all_lineages.append(top_lineage)
                    seen_lineages.add(top_lineage)

            for s_name, sample_data in samples.items():
                for lineage, _ in sorted(sample_data.items(), key=lambda xv: xv[1], reverse=True):
                    if lineage not in seen_lineages:
                        all_lineages.append(lineage)
                        seen_lineages.add(lineage)

            self.general_stats_cols(top_lineages_dict, all_lineages, taxid, organism)
            self.add_freyja_section(all_lineages, samples, taxid, organism)

        # Write combined data file at the very end
        combined: Dict[str, Dict] = {}
        for taxid, samples in data_by_taxid.items():
            organism = PATHOGEN_MAP[taxid]
            for s_name, data in samples.items():
                combined[f"{s_name} ({organism})"] = data
        self.write_data_file(combined, "multiqc_freyja")

    def general_stats_cols(self, top_lineages_dict, all_lineages, taxid, organism):
        """Add a pair of columns per organism to the General Statistics table."""
        headers = {
            f"Top_lineage_freyja_{taxid}": {
                "title": f"Top lineage ({organism})",
                "description": f"The most abundant lineage in the sample ({organism})",
                "bgcols": {x: self.scale.get_colour(i) for i, x in enumerate(all_lineages)},
            },
            f"Top_lineage_freyja_{taxid}_percentage": {
                "title": f"Top lineage % ({organism})",
                "description": f"The percentage of the most abundant lineage in the sample ({organism})",
                "max": 100,
                "min": 0,
                "scale": "Blues",
                "modify": lambda x: x * 100,
                "suffix": "%",
            },
        }
        self.general_stats_addcols(top_lineages_dict, headers)

    def add_freyja_section(self, lineages, data_by_sample, taxid, organism):
        pconfig = {
            "id": f"Freyja_plot_{taxid}",
            "title": f"Freyja: Top lineages ({organism})",
            "ylab": "relative abundance",
            "y_clipmax": 1,
            "cpswitch": False,
            "cpswitch_c_active": False,
        }
        cats = {x: {"name": x, "color": self.scale.get_colour(i, lighten=1)} for i, x in enumerate(lineages)}

        self.add_section(
            name=f"Freyja Summary ({organism})",
            anchor=f"freyja-summary-{taxid}",
            description="""
                Relative lineage abundances from mixed samples. Hover over the column headers for descriptions and click _Help_ for more in-depth documentation.
                """,
            helptext="""
                The graph denotes a sum of all lineage abundances in a particular WHO designation, otherwise they are grouped into "Other".
                Lineages abundances are calculated as the number of reads that are assigned to a particular lineage.
                Lineages and their corresponding abundances are summarized by constellation.

                > **Note**: Lineage designation is based on the used WHO nomenclature, which could vary over time.
                """,
            plot=bargraph.plot(data_by_sample, cats, pconfig),
        )
