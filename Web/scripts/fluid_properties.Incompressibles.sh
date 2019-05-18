#!/bin/bash
set -e
SEP="####################################"
echo $SEP
echo "Running a nightly build on $(date)"
echo $SEP
pushd ../../dev/incompressible_liquids
echo "Generating fitting reports"
python all_incompressibles.py -ns
echo "Creating example figures"
pdfnup --quiet --nup 2x1 --delta "1cm 0cm" -o report/report2up.pdf report/DowQ_fitreport.pdf report/LiBr_fitreport.pdf
convert -background "#FFFFFF" -density 300 report/report2up.pdf report/report2up.jpg # Convert the PDF to JPG
convert -crop 100%x47%+0+30 -resize '850x' -quality 75 report/report2up.jpg report/report2up.jpg # Resize it
echo "Copying the reports to Web/_static/fluid_properties"
mkdir -p ../../Web/_static/fluid_properties/Incompressibles_reports
rsync -am report/ ../../Web/_static/fluid_properties/Incompressibles_reports
echo "Copying the tables to Web/fluid_properties"
mkdir -p ../../Web/fluid_properties
rsync -am tables/ ../../Web/fluid_properties
popd
echo $SEP
echo "Generating validation figures"
python fluid_properties.Incompressibles.py
echo "Resizing figures and creating JPG"
convert -background "#FFFFFF" -density 72 incompressibles_consistency.pdf incompressibles_consistency.jpg # Convert the PDF to JPG
mkdir -p ../../Web/_static/fluid_properties
cp incompressibles_consistency.pdf incompressibles_consistency.jpg ../../Web/_static/fluid_properties
echo $SEP
echo "Done"
exit 0
