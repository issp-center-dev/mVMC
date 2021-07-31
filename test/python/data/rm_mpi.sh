for name in HeisenbergChain HeisenbergChain_cmp HeisenbergChain_fsz \
HubbardChain HubbardChain_cmp HubbardChain_fsz HubbardTetragonal HubbardTetragonal_MomentumProjection \
KondoChain KondoChain_cmp KondoChain_fsz KondoChain_Stot1_cmp
do
echo "delete $name\_mpi"
rm -r $name\_mpi
done
