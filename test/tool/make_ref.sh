x=1
end=11
rm -r ./output*
while [ $x -lt $end ]
do
    x0=`expr 5 \* $x + 12345`
    cp ./StdFace_org.def ./StdFace.def
    echo "RndSeed = $x0" >> StdFace.def
    ./vmc.out -s ./StdFace.def ./initial.def
    mv ./output ./output$x
    x=`expr $x + 1`
done

cp ./StdFace_org.def ./StdFace_org.def.bak
rm ./*.def
rm ./*.dat
rm ./*.gp
mv ./StdFace_org.def.bak ./StdFace_org.def

