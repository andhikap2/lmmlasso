#EMMAX


./plink --bfile k2029 --recode transpose --output-missing-genotype 0 --out k2029_emmax 

./emmax-kin-intel64 -v -s -d 10 k2029_emmax  #IBS matrix

./emmax-kin-intel64 -v -d 10 k2029_emmax #Balding-Nichols Matrix 


#EMMAX k2030 matrix: add an extra line to the tfam file, then read the tped, add an extra row of randomly generated 0's and 1's, and see how the overall scaling changes... Derive some sort of linear solver of that...

#Stepd

'''

1) add a new observation with a randomly generated SNP dataset
2) recompute the EMMAX IBS
3) subset the original k2029 from k3030, call it k2029*
4) figure out the linear relationship between k2029 and k2029*


notes:
reading tped file requires installing pandas-plink

from os.path import join
from pandas_plink import read_plink
from pandas_plink import get_data_folder
tped=read_plink(join(get_data_folder(),"k2029_emmax"))
'''



./plink --bfile k2030_test --recode 12 transpose --output-missing-genotype 0 --out k2030_emmax 

./emmax-kin-intel64 -v -s -d 10 k2030_emmax  #IBS matrix

./emmax-kin-intel64 -v -d 10 k2029_emmax #Balding-Nichols Matrix 
