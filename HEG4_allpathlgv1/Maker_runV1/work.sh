python ../../Maker.py --config test.config
perl /rhome/cjinfeng/BigData/00.RD/Annotation/HEG4/Maker/bin/qsub-pbs.pl --workdir `pwd`/Maker.run.cut --convert no --lines 2 --queue js --resource nodes=1:ppn=1,walltime=200:00:00 maker.sh
perl /rhome/cjinfeng/BigData/00.RD/Annotation/HEG4/Maker/bin/qsub-pbs.pl --workdir `pwd`/Maker.run.cut --convert no maker.gff3.sh

