1) Compilation:
cd  work_2.2.hybrid
ls *ake*
makefile.haswell  makefile.knl
ex: 
make -f makefile.knl clean
rm *.o  *.mod
make -f makefile.knl
mpif90  -extend-source 132  -xMIC-AVX512 -O3 -c module_don.f
mpif90   -extend-source 132 -xMIC-AVX512 -O3 -c module_subr.f
mpif90  -extend-source 132 -xMIC-AVX512 -O3 -c prog_pn.f
mpif90  -DOPENMP -extend-source 132 -xMIC-AVX512 -O3 -c CALCUL.F90  -qopenmp
mpif90   -extend-source 132 -xMIC-AVX512 -O3 -o Xpn_knl module_don.o
module_subr.o  prog_pn.o CALCUL.o          -qopenmp

2) Run

cd RUN_BIG_LONG (800s on knl in sequential mode, 200s on haswell)
cat sub_knl , adapt to your system launcher
ccc_msub sub_knl

numerical result are written in fort.440 

3) Check result if necessary:
For checking result availability, code diagnostic compares fort.440 to the one 
of REF/res_BIG_LONG.tar.gz
cd REF
tar xvfz REF/res_BIG_LONG.tar.gz 
mv fort.440 fort.410 (if no file fort.410)
cp fort.410 ../RUN_BIG_LONG
cd ../RUN_BIG_LONG
./diagnostic
Je lis fort.440 matrice de test
Je lis fort.410 matrice de reference
   DIFFERENCE =   7.264978549303969E-014
   somme de toutes les diff = -1.451986077160262E-013

4) INFO

 The hotspot of the miniapp is 
function RESC_pn(valeur12)  in CALCUL.F90
This fuacntion will be called for each matrix element 
and the number of call is preprocessor dependant.
In real life, the number of calls will be very high.

For test/modif purpose we have reduce the run to 800s 
(80s of init) on KNL in RUN_BIG_LONG 

function prepa_vabcd will be called only one time 
at the begin so can be ignored





