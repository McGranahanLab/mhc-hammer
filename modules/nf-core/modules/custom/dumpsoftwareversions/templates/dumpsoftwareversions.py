#!/usr/bin/env python

# This is a modified version of the original dumpsoftwareversions.py nf-core script


"""Provide functions to merge multiple versions.yml files."""


import yaml
import platform

def main():
    """Load all version files and generate merged output."""
    versions_this_module = {}
    versions_this_module["NF_MHC_HAMMER:MHC_HAMMER:CUSTOM_DUMPSOFTWAREVERSIONS"] = {
        "python": platform.python_version(),
        "yaml": yaml.__version__,
    }

    with open("collated_versions.yml") as f:
        versions_by_process = yaml.load(f, Loader=yaml.BaseLoader) | versions_this_module

    # aggregate versions by the module name (derived from fully-qualified process name)
    versions_by_module = {}
    for process, process_versions in versions_by_process.items():
        module = process.split(":")[-1]
        try:
            if versions_by_module[module] != process_versions:
                raise AssertionError(
                    "We assume that software versions are the same between all modules. "
                    "If you see this error-message it means you discovered an edge-case "
                    "and should open an issue. "
                )
        except KeyError:
            versions_by_module[module] = process_versions

    versions_by_module["Workflow"] = {
        "Nextflow": "21.10.6",
        "mhc_hammer": "1.0",
    }

    with open("software_versions.yml", "w") as f:
        yaml.dump(versions_by_module, f, default_flow_style=False)
    with open("versions.yml", "w") as f:
        yaml.dump(versions_this_module, f, default_flow_style=False)


if __name__ == "__main__":
    main()