# SLAT - <ins>S</ins>hotgun <ins>L</ins>ipidomic <ins>A</ins>ssignment <ins>T</ins>ool #


[This webapp](https://briankleiboeker.shinyapps.io/structure_from_comp/) was designed for use with high-resolution Orbitrap ESI-MS data with the goal of accurately assigining lipid class and structure to user-provided elemental compositions or m/z values, but should also be compatible with most any mass spectromerty data so long as the elemental compositions and/or column names are provided in a standard format. 

For high-throughput, offline, or comparative analyses, check out the associated R package, [slatR](https://github.com/briankleiboeker/slatR).
![abstract](https://github.com/briankleiboeker/SLAT/assets/59810795/44ae5690-83a4-4fce-b2ef-f977d0e358bb)

To use SLAT to assign lipid structure to a single elemental composition or m/z value, see this youtube video:

[![singlecompvid](http://img.youtube.com/vi/FBKgMt7WbcY/0.jpg)](https://youtu.be/FBKgMt7WbcY "Video Title")

For muliple elemental compositions or m/z values (a test dataset is available in this repository):

[![multicompvid](http://img.youtube.com/vi/Bjx4WziZP0c/0.jpg)](https://youtu.be/Bjx4WziZP0c "Video Title")

To use SLAT to recalibrate raw data (from m/z values):

[![recalibrationvid](http://img.youtube.com/vi/QjfLFrz4-TA/0.jpg)](https://youtu.be/QjfLFrz4-TA "Video Title")


## Change log ##
* v1.0.0 Release SLAT
* v1.0.1 Add exact structure database. SLAT now returns most likely exact structure(s) for over 200 general structures. The exact structure database is based on hundreds of MS2 and MS3 analyses of lipids from various tissues and cell lines.
* v1.0.2 Add ability to assign elemental composition to m/z values and add 3 specific bacterial strains to exact structure database (e.coli, streptococcus, and listeria).
