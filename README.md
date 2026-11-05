<p align="center">
  <img src="inst/www/mamamia-logo.svg" alt="MaMaMIA logo" width="96" height="96">
</p>

<h1 align="center">MaMaMIA</h1>

<p align="center">
  <strong>Max Mapping for Meticulous Introgression Analysis</strong><br>
  Version 0.1.0 · GPL (>= 3)
</p>

MaMaMIA is an interactive dashboard for detecting interspecific introgressions
from Genotyping-by-Sequencing (GBS) coverage. Upload a coordinate-sorted BAM,
choose a donor and recipient species, run the analysis, and explore the results
as an interactive karyotype, sortable call table, and downloadable report
bundle. The method is ploidy-agnostic and supports any number of donor and
recipient components.

## Quick start

1. Open the MaMaMIA dashboard provided by your deployment.
2. Upload a BAM, or choose **Example data** for a guided tour.
3. Select the recipient and donor species.
4. Pick an analysis profile and press **Run analysis**.
5. Explore **Karyotype** and **Calls**, then download the result bundle from
   **Export**.

For BAMs above the deployment upload limit, slim them locally with
[bamhusk](https://github.com/alermol/bamhusk) before uploading.

## Documentation

Detailed guides, screenshots, configuration, method, and deployment notes live
in the [MaMaMIA wiki](https://github.com/alermol/MaMaMIA/wiki).

| Guide | Contents |
| --- | --- |
| [Getting Started](https://github.com/alermol/MaMaMIA/wiki/Getting-Started) | First run, Docker, and source setup |
| [Dashboard Guide](https://github.com/alermol/MaMaMIA/wiki/Dashboard-Guide) | Every control, filter, tab, and download |
| [Inputs and Assemblies](https://github.com/alermol/MaMaMIA/wiki/Inputs-and-Assemblies) | BAM requirements, bamhusk, registry, components, and pairing |
| [Analysis Profiles](https://github.com/alermol/MaMaMIA/wiki/Analysis-Profiles) | Balanced, Sensitive, CNV robust, Complex, and Custom |
| [Method](https://github.com/alermol/MaMaMIA/wiki/Method) | Windows, scoring, masks, CBS, calibration, and placement |
| [Results and Calls](https://github.com/alermol/MaMaMIA/wiki/Results-and-Calls) | Reading classes, status, placement, and dosage |
| [Export Bundle](https://github.com/alermol/MaMaMIA/wiki/Export-Bundle) | Files, formats, QC, and reproducibility |
| [Configuration Reference](https://github.com/alermol/MaMaMIA/wiki/Configuration-Reference) | Complete registry schema and defaults |
| [Deployment](https://github.com/alermol/MaMaMIA/wiki/Deployment) | Docker, ShinyProxy, mounts, limits, and security |
| [Troubleshooting](https://github.com/alermol/MaMaMIA/wiki/Troubleshooting) | Common failures and fixes |

## License

GPL (>= 3). See [`LICENSE.md`](LICENSE.md).

## Contact

Aleksey Ermolaev <alerm@al-erm.ru> ·
[Issues](https://github.com/alermol/MaMaMIA/issues)
