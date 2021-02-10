!/bin/bash
#

for SEASONALITY in `seq 1 2`
do 
for DECAY in `seq 1 3`
do
for ACESS in `seq 1 2`
do
for AGE in 1 
do
for EIR in `seq 1 1 5`
do
    sbatch launch_gp.sh $SEASONALITY $DECAY $ACESS $AGE $EIR
done
done
done
done
done
