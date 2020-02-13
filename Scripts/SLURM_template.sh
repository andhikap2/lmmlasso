#!/bin/bash
#Making microclimate files

k=0
for TEMP in 0  ; do for PHOTO in 0 ; do for PERIOD in 24; do for RSQUARE in 0.1 K ; do for KFILE in *K_MATRIX**.npy*; do k=$(($k+1))
sed "s/RSQUARE/$RSQUARE/" dtblmmlasso_testk.slurm > dtb$k.slurm
sed -i "s/TEMP/$TEMP/" dtb$k.slurm
sed -i "s/PERIOD/$PERIOD/" dtb$k.slurm
sed -i "s/PHOTO/$PHOTO/" dtb$k.slurm
sed -i "s/KFILE/$KFILE/" dtb$k.slurm
done
done
done
done
done

for i in *[0-9].slurm; do sbatch $i
done