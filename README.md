# SLAT - <ins>S</ins>hotgun <ins>L</ins>ipidomic <ins>A</ins>ssignment <ins>T</ins>ool #


[This webapp](https://lodhilab.shinyapps.io/slat/) was designed for use with high-resolution Orbitrap ESI-MS data with the goal of accurately assigining lipid class and structure to user-provided elemental compositions, but should also be compatible with most any mass spectromerty data so long as the elemental compositions and column names are provided in a standard format. 

For high-throughput, offline, or comparative analyses, check out the associated R package, [slatR](https://github.com/briankleiboeker/slatR).
![abstract](https://github.com/briankleiboeker/SLAT/assets/59810795/44ae5690-83a4-4fce-b2ef-f977d0e358bb)
To use SLAT to assign lipid structure to a single elemental composition:

![single composition](single_comp_demo.gif)

For muliple elemental compositions (a test dataset is available in this repository):

![multiple compositions](bulk_comps_demo.gif)

## Change log ##
* v1.0.0 Release SLAT
* v1.0.1 Add exact structure database. SLAT now returns most likely exact structure(s) for over 200 general structures. The exact structure database is based on hundreds of MS2 and MS3 analyses of lipids from various tissues and cell lines.
