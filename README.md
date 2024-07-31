# BroadBandSeabed
Broad-band frequency analysis of underwater acoustic seabed response

## Organization of the repository
* `/RAW/`: Measurements obtained from SimRad EK80.
* `/FIG/`: Figures generated for the thesis report.
* `/ESP3/`: Code obtained from ESP3 for generating Matlab compatible data.
* `/gps_position.m`: Derive the GPS positions of every ping.
* `/example.m`: Companion code for Section 3 of the report.
* `/precomputation_absorption.m`: Precomputation for the absorption losses.
* `/precomputation_pulse.m`: Precomputation for getting the input pulse.
* `/impulse_response.m`: Compute the corrections on the signals and the average impulse response for the different pings.
* `/bbc.m`: File for the characterization of the sediment.
* `/bbd.m`: File for the detection of buried objects.

## Data from the measurements
The measurements are split into files of around 150Mb size. We do not provide the 158 files, but you can ask them by email at `simon.masson@protonmail.com`. One file is provided for running `example.m`.

