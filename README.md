# WOCE_GLODAP_Conversions
Conversions from WOCE Section Names to expocodes used in GOSHIP and GLODAP

Each $SEC.txt file contains nominal years for an occupation, as well as expocodes
for both GO-SHIP and GLODAP. GO-SHIP expocodes have been obtained from [CCDHO](https://cchdo.ucsd.edu/)
and GLODAP expocodes obtained from [NCEI](https://www.ncei.noaa.gov/access/ocean-carbon-data-system/oceans/RepeatSections/).
These codes occasionally differ: this project is an attempt to categorise these
differences.

This list is a work in progress and will continue to grow. The idea is to create
a compatibility layer which allows the use of WOCE Hydrographic sections in 
conjuction with tools like [GO-SHIP-Easy-Ocean](https://github.com/kkats/GO-SHIP-Easy-Ocean) to easily access biogeochemical 
data from GLODAP for repeat occupations, without having to use tools such as ODV.

At present, data should be treated with caution, as these lists of occupations for
each hydrographic section have not been checked, nor has any extra data processing
been performed to account for partial occupations or cruises which took data on 
more than one WOCE hydrographic section. These issues will be resolved over time.
For now these data should be treated as a starting point rather than a canonical
guide.
