python processor/run.py processor/cfg/xcone.cfg
rm -r processor/histos
rm -r plots/
mkdir plots
mv histos/ processor/histos
cd processor/
python simple_plot.py
cd ../
