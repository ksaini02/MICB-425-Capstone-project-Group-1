mamba activate treesapp_cenv
treesapp assign -i MICB425/assemblies/SI072_10m_contig.fa -o ~/TS_NapA_10m_output_[group]_[init]/ -t NapA -n 2
treesapp assign -i MICB425/assemblies/SI072_100m_contig.fa -o ~/TS_NapA_100m_output_[group]_[init]/ -t NapA -n 2
treesapp assign -i MICB425/assemblies/SI072_120m_contig.fa -o ~/TS_NapA_120m_output_[group]_[init]/ -t NapA -n 2
treesapp assign -i MICB425/assemblies/SI072_135m_contig.fa -o ~/TS_NapA_135m_output_[group]_[init]/ -t NapA -n 2
treesapp assign -i MICB425/assemblies/SI072_150m_contig.fa -o ~/TS_NapA_150m_output_[group]_[init]/ -t NapA -n 2
treesapp assign -i MICB425/assemblies/SI072_165m_contig.fa -o ~/TS_NapA_165m_output_[group]_[init]/ -t NapA -n 2
treesapp assign -i MICB425/assemblies/SI072_200m_contig.fa -o ~/TS_NapA_200m_output_[group]_[init]/ -t NapA -n 2

guppy fpd -o [TS_output]/iTOL_output/NapA/alpha_diversiy.csv --csv [TS_output]/iTOL_output/NapA/NapA_complete_profile.jplace 