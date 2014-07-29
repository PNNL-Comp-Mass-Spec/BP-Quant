clear

cd \\PNL\Projects\ProteomicsCollab\iTRAQ
protein_list = {'sll1579' 'sll1580', 'sll1578', 'sll1577', 'slr2067', 'slr1986', 'slr0335'};

load (protein_list{i})
d = nan(size(ds));
d(:,4:15) = double(ds(:,4:15));
CB_ratio = [d(:,4)./d(:,7), d(:,8)./d(:,11), d(:,12)./d(:,15)];
CK_ratio = [d(:,5)./d(:,7), d(:,9)./d(:,11), d(:,13)./d(:,15)];
PAL_ratio = [d(:,6)./d(:,7), d(:,10)./d(:,11), d(:,14)./d(:,15)];
PROTEIN_SIG = [nanmean(CB_ratio,2), nanmean(CK_ratio,2), nanmean(PAL_ratio,2)];
    
cd '\\PNL\Projects\ProteomicToolbox\BPQuant\BP_Quant Toolbox'
[ PEPTIDE_IDX, POST_PROB, NUM_PROTEOFORMS] = BPQuant( PROTEIN_SIG, PI_NOT );