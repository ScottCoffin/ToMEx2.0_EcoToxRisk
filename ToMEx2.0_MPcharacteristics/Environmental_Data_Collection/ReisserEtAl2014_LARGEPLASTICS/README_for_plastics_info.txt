====================
readme file of Reisser et al (2013) Data from: Marine plastic pollution in waters around Australia: characteristics, concentrations, and pathways. Dryad Digital Repository. 

doi:10.5061/dryad.mt0j5

Feedback: jureisser@gmail.com

================================
Description of files

readme.txt --> this document
nettows_info.csv --> contains information related to the sea surface net tows undertaken (columns = variables, lines = net tows)
plastics_info.csv --> contains information related to the marine plastics collected during the study (columns = variables, lines = plastic pieces)

================================
Description of the variables of nettows_info.csv matrix (see its heading/first line)

VesselTripName --> vessel and voyage name
NetStation --> net station identification number (where 3 consecutive net tows were undertaken); goes from 1 to 57
Replicate --> indicates if it is the first (1), second (2), or third (3) net tow undertaken at the net station
CollectedBy --> states who wrote down the variables of nettows_info.csv and plastics_info.csv in our field datasheets 
DateUTC --> date in coordinated universal time; day.month.year
StartTimeUTC --> starting time in coordinated universal time; hour:minute
EndTimeUTC --> ending time in coordinated universal time; hour:minute
Duration --> duration in minutes of each net tow
NetTypeMeshSizeDimensionsNetMouth --> indicate net type (Manta or Neuston), mesh size (in microns), and net frame dimensions (width x height, in m)
VesselSpeed --> vessel speed (in knots)
VesselDirection --> vessel direction (in degrees)
WindSpeed --> wind speed (in knots)
WindDirection --> wind direction (in degrees)
SeaTemperature --> sea surface temperature (in degrees Celsius)
Salinity --> sea surface salinity
Depth --> local depth (in meters)
StartLatitute --> starting Latitute; in degrees
StartLongitude --> starting Longitude; in degrees
EndLatitude --> ending Latitute; in degrees
EndLongitude --> ending Longitude; in degrees
StartFlowmeter --> starting flowmeter value (see flowmeter details here: http://envco.com.au/catalog/product/open-channel-flow/general-oceanics-mechanical-flowmeter.html)
EndFlowmeter --> ending flowmeter value
Flow --> total flow (ending flowmeter value minus starting flowmeter value)
HardPlastics --> number of hard plastic fragments
SoftPlastics --> number of soft plastic fragments
PlasticLines --> number of plastic line fragments
Styrofoam --> number of styrofoam fragments
Pellets --> number of plastic pellets/nurdles
TotalPlastics --> total number of plastic pieces (all types)
DriftingWood --> number of floating wood fragments
Pumice --> number of pumice particles
Cs --> sea surface plastic concentration (pieces per km^2)

================================
Description of the variables of plastics_info.csv matrix (see its heading/first line)

Plastic --> plastic ID number
NetStation --> net station identification number (where 3 consecutive net tows were undertaken); goes from 1 to 57
Replicate --> indicate if it is from the first (1), second (2), or third (3) net tow undertaken at the net station
Type --> type of plastic (e.g. hard, soft, styrofoam)
Length --> plastic size (in mm)
Colour --> plastic colour
PolymerType --> plastic polymer as determined by FT-IR spectra


============
Acknowledgements

We thank (1) CSIRO Marine National Facility and AIMS data centre for the 'underway' (environmental) data; (3) the International Postgraduate Research Scholarship (IPRS)  and the CSIRO Postgraduate scholarship program for funding; (3) CSIRO, AIMS, and Austral Fisheries for seatime aboard their vessels.

==================
end of readme file
