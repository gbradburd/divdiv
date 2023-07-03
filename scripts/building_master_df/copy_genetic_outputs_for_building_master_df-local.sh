#!/bin/bash

echo starting dataset bioprj_PRJNA371502_Acanthurus-olivaceus
rclone copy divdiv_datafiles:bioprj_PRJNA553831_Acanthephyra-purpurea/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA371502_Acanthurus-olivaceus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA473816_Acropora-cervicornis
rclone copy divdiv_datafiles:bioprj_PRJNA371502_Acanthurus-olivaceus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA473816_Acropora-cervicornis/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA593014_Acropora-millepora
rclone copy divdiv_datafiles:bioprj_PRJNA473816_Acropora-cervicornis/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA593014_Acropora-millepora/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA473816_Acropora-palmata
rclone copy divdiv_datafiles:bioprj_PRJNA593014_Acropora-millepora/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA473816_Acropora-palmata/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA473816_Acropora-prolifera
rclone copy divdiv_datafiles:bioprj_PRJNA473816_Acropora-palmata/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA473816_Acropora-prolifera/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA361144_Agaricia-fragilis
rclone copy divdiv_datafiles:bioprj_PRJNA473816_Acropora-prolifera/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA361144_Agaricia-fragilis/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA385083_Agaricia-undata
rclone copy divdiv_datafiles:bioprj_PRJNA361144_Agaricia-fragilis/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA385083_Agaricia-undata/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA294760_Amphiprion-bicinctus
rclone copy divdiv_datafiles:bioprj_PRJNA385083_Agaricia-undata/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA294760_Amphiprion-bicinctus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA563695_Amphiprion-clarkii
rclone copy divdiv_datafiles:bioprj_PRJNA294760_Amphiprion-bicinctus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA563695_Amphiprion-clarkii/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA295681_Anguilla-anguilla
rclone copy divdiv_datafiles:bioprj_PRJNA563695_Amphiprion-clarkii/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA295681_Anguilla-anguilla/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA295681_Anguilla-rostrata
rclone copy divdiv_datafiles:bioprj_PRJNA295681_Anguilla-anguilla/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA295681_Anguilla-rostrata/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA632874_Anthopleura-elegantissima
rclone copy divdiv_datafiles:bioprj_PRJNA295681_Anguilla-rostrata/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA632874_Anthopleura-elegantissima/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA436919_Apostichopus-californicus
rclone copy divdiv_datafiles:bioprj_PRJNA632874_Anthopleura-elegantissima/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA436919_Apostichopus-californicus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA343012_Aptenodytes-patagonicus
rclone copy divdiv_datafiles:bioprj_PRJNA436919_Apostichopus-californicus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA343012_Aptenodytes-patagonicus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA419445_Arctocephalus-forsteri
rclone copy divdiv_datafiles:bioprj_PRJNA343012_Aptenodytes-patagonicus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA419445_Arctocephalus-forsteri/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA299372_Aythya-marila
rclone copy divdiv_datafiles:bioprj_PRJNA419445_Arctocephalus-forsteri/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA299372_Aythya-marila/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA542967_Bartholomea-annulata
rclone copy divdiv_datafiles:bioprj_PRJNA299372_Aythya-marila/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA542967_Bartholomea-annulata/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA595396_Bathygobius-cocosensis
rclone copy divdiv_datafiles:bioprj_PRJNA542967_Bartholomea-annulata/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA595396_Bathygobius-cocosensis/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA473646_Bathymodiolus-platifrons
rclone copy divdiv_datafiles:bioprj_PRJNA595396_Bathygobius-cocosensis/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA473646_Bathymodiolus-platifrons/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA646172_Bathyraja-aleutica
rclone copy divdiv_datafiles:bioprj_PRJNA473646_Bathymodiolus-platifrons/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA646172_Bathyraja-aleutica/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA646172_Bathyraja-panthera
rclone copy divdiv_datafiles:bioprj_PRJNA646172_Bathyraja-aleutica/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA646172_Bathyraja-panthera/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA646172_Bathyraja-parmifera
rclone copy divdiv_datafiles:bioprj_PRJNA646172_Bathyraja-panthera/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA646172_Bathyraja-parmifera/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA340157_Callinectes-sapidus
rclone copy divdiv_datafiles:bioprj_PRJNA646172_Bathyraja-parmifera/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA340157_Callinectes-sapidus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA559677_Caretta-caretta
rclone copy divdiv_datafiles:bioprj_PRJNA340157_Callinectes-sapidus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA559677_Caretta-caretta/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA286949_Chlorostoma-funebralis
rclone copy divdiv_datafiles:bioprj_PRJNA559677_Caretta-caretta/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA286949_Chlorostoma-funebralis/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA430030_Coryphaenoides-rupestris
rclone copy divdiv_datafiles:bioprj_PRJNA286949_Chlorostoma-funebralis/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA430030_Coryphaenoides-rupestris/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA573759_Cranchia-scabra
rclone copy divdiv_datafiles:bioprj_PRJNA430030_Coryphaenoides-rupestris/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA573759_Cranchia-scabra/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA508986_Ctenolabrus-rupestris
rclone copy divdiv_datafiles:bioprj_PRJNA573759_Cranchia-scabra/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA508986_Ctenolabrus-rupestris/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA590458_Dascyllus-trimaculatus
rclone copy divdiv_datafiles:bioprj_PRJNA508986_Ctenolabrus-rupestris/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA590458_Dascyllus-trimaculatus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA564823_Delphinapterus-leucas
rclone copy divdiv_datafiles:bioprj_PRJNA590458_Dascyllus-trimaculatus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA564823_Delphinapterus-leucas/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA420851_Ectopleura-larynx
rclone copy divdiv_datafiles:bioprj_PRJNA564823_Delphinapterus-leucas/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA420851_Ectopleura-larynx/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA311981_Engraulis-encrasicolus
rclone copy divdiv_datafiles:bioprj_PRJNA420851_Ectopleura-larynx/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA311981_Engraulis-encrasicolus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA523574_Eudyptes-chrysolophus
rclone copy divdiv_datafiles:bioprj_PRJNA311981_Engraulis-encrasicolus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA523574_Eudyptes-chrysolophus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA576675_Eukrohnia-hamata
rclone copy divdiv_datafiles:bioprj_PRJNA523574_Eudyptes-chrysolophus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA576675_Eukrohnia-hamata/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA394157_Exaiptasia-diaphana
rclone copy divdiv_datafiles:bioprj_PRJNA576675_Eukrohnia-hamata/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA394157_Exaiptasia-diaphana/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJDB7819_Fibramia-amboinensis
rclone copy divdiv_datafiles:bioprj_PRJNA394157_Exaiptasia-diaphana/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJDB7819_Fibramia-amboinensis/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA629489_Fucus-vesiculosus
rclone copy divdiv_datafiles:bioprj_PRJDB7819_Fibramia-amboinensis/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA629489_Fucus-vesiculosus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA437462_Fundulus-grandis
rclone copy divdiv_datafiles:bioprj_PRJNA629489_Fucus-vesiculosus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA437462_Fundulus-grandis/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA477712_Fundulus-heteroclitus
rclone copy divdiv_datafiles:bioprj_PRJNA437462_Fundulus-grandis/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA477712_Fundulus-heteroclitus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA558810_Gadus-macrocephalus
rclone copy divdiv_datafiles:bioprj_PRJNA477712_Fundulus-heteroclitus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA558810_Gadus-macrocephalus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA528403_Gadus-morhua
rclone copy divdiv_datafiles:bioprj_PRJNA558810_Gadus-macrocephalus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA528403_Gadus-morhua/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA576132_Galaxea-fascicularis
rclone copy divdiv_datafiles:bioprj_PRJNA528403_Gadus-morhua/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA576132_Galaxea-fascicularis/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA576132_Galaxea-horrescens
rclone copy divdiv_datafiles:bioprj_PRJNA576132_Galaxea-fascicularis/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA576132_Galaxea-horrescens/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA564121_Galaxias-maculatus
rclone copy divdiv_datafiles:bioprj_PRJNA576132_Galaxea-horrescens/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA564121_Galaxias-maculatus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA453151_Gasterosteus-aculeatus
rclone copy divdiv_datafiles:bioprj_PRJNA564121_Galaxias-maculatus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA453151_Gasterosteus-aculeatus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA448438_Holacanthus-passer
rclone copy divdiv_datafiles:bioprj_PRJNA453151_Gasterosteus-aculeatus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA448438_Holacanthus-passer/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA281764_Homarus-americanus
rclone copy divdiv_datafiles:bioprj_PRJNA448438_Holacanthus-passer/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA281764_Homarus-americanus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA551577_Idotea-baltica
rclone copy divdiv_datafiles:bioprj_PRJNA281764_Homarus-americanus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA551577_Idotea-baltica/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA643849_Isocladus-armatus
rclone copy divdiv_datafiles:bioprj_PRJNA551577_Idotea-baltica/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA643849_Isocladus-armatus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA553554_Istiblennius-lineatus
rclone copy divdiv_datafiles:bioprj_PRJNA643849_Isocladus-armatus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA553554_Istiblennius-lineatus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA280898_Lagenorhynchus-acutus
rclone copy divdiv_datafiles:bioprj_PRJNA553554_Istiblennius-lineatus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA280898_Lagenorhynchus-acutus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA280898_Lagenorhynchus-albirostris
rclone copy divdiv_datafiles:bioprj_PRJNA280898_Lagenorhynchus-acutus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA280898_Lagenorhynchus-albirostris/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA480409_Laguncularia-racemosa
rclone copy divdiv_datafiles:bioprj_PRJNA280898_Lagenorhynchus-albirostris/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA480409_Laguncularia-racemosa/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA533707_Larimichthys-polyactis
rclone copy divdiv_datafiles:bioprj_PRJNA480409_Laguncularia-racemosa/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA533707_Larimichthys-polyactis/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA356786_Lateolabrax-maculatus
rclone copy divdiv_datafiles:bioprj_PRJNA533707_Larimichthys-polyactis/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA356786_Lateolabrax-maculatus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA521302_Leptopsammia-pruvoti
rclone copy divdiv_datafiles:bioprj_PRJNA356786_Lateolabrax-maculatus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA521302_Leptopsammia-pruvoti/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA329407_Lutjanus-campechanus
rclone copy divdiv_datafiles:bioprj_PRJNA521302_Leptopsammia-pruvoti/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA329407_Lutjanus-campechanus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJDB7819_Lutjanus-fulvus
rclone copy divdiv_datafiles:bioprj_PRJNA329407_Lutjanus-campechanus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJDB7819_Lutjanus-fulvus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA631144_Mallotus-villosus
rclone copy divdiv_datafiles:bioprj_PRJDB7819_Lutjanus-fulvus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA631144_Mallotus-villosus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA544861_Montipora-capitata
rclone copy divdiv_datafiles:bioprj_PRJNA631144_Mallotus-villosus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA544861_Montipora-capitata/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA564770_Mytilus-galloprovincialis
rclone copy divdiv_datafiles:bioprj_PRJNA544861_Montipora-capitata/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA564770_Mytilus-galloprovincialis/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA369717_Nautilus-pompilius
rclone copy divdiv_datafiles:bioprj_PRJNA564770_Mytilus-galloprovincialis/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA369717_Nautilus-pompilius/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA430897_Nymphon-australe
rclone copy divdiv_datafiles:bioprj_PRJNA369717_Nautilus-pompilius/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA430897_Nymphon-australe/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA343301_Oculina-patagonica
rclone copy divdiv_datafiles:bioprj_PRJNA430897_Nymphon-australe/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA343301_Oculina-patagonica/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA508589_Orbicella-faveolata
rclone copy divdiv_datafiles:bioprj_PRJNA343301_Oculina-patagonica/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA508589_Orbicella-faveolata/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA511386_Ostrea-lurida
rclone copy divdiv_datafiles:bioprj_PRJNA508589_Orbicella-faveolata/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA511386_Ostrea-lurida/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA373825_Paracentrotus-lividus
rclone copy divdiv_datafiles:bioprj_PRJNA511386_Ostrea-lurida/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA373825_Paracentrotus-lividus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA640135_Paracirrhites-forsteri
rclone copy divdiv_datafiles:bioprj_PRJNA373825_Paracentrotus-lividus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA640135_Paracirrhites-forsteri/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA477007_Paralichthys-dentatus
rclone copy divdiv_datafiles:bioprj_PRJNA640135_Paracirrhites-forsteri/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA477007_Paralichthys-dentatus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA554161_Penaeus-duorarum
rclone copy divdiv_datafiles:bioprj_PRJNA477007_Paralichthys-dentatus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA554161_Penaeus-duorarum/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA450328_Halichoerus-grypus-atlantica
rclone copy divdiv_datafiles:bioprj_PRJNA554161_Penaeus-duorarum/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA450328_Halichoerus-grypus-atlantica/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA343959_Phocoena-phocoena
rclone copy divdiv_datafiles:bioprj_PRJNA450328_Halichoerus-grypus-atlantica/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA343959_Phocoena-phocoena/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA659918_Phocoena-sinus
rclone copy divdiv_datafiles:bioprj_PRJNA343959_Phocoena-phocoena/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA659918_Phocoena-sinus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA445895_Pisaster-ochraceus
rclone copy divdiv_datafiles:bioprj_PRJNA659918_Phocoena-sinus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA445895_Pisaster-ochraceus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA471559_Planes-minutus
rclone copy divdiv_datafiles:bioprj_PRJNA445895_Pisaster-ochraceus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA471559_Planes-minutus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA382467_Platichthys-flesus
rclone copy divdiv_datafiles:bioprj_PRJNA471559_Planes-minutus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA382467_Platichthys-flesus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA450871_Pocillopora-damicornis
rclone copy divdiv_datafiles:bioprj_PRJNA382467_Platichthys-flesus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA450871_Pocillopora-damicornis/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA280774_Poecilia-latipinna
rclone copy divdiv_datafiles:bioprj_PRJNA450871_Pocillopora-damicornis/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA280774_Poecilia-latipinna/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA559876_Pomacanthus-maculosus
rclone copy divdiv_datafiles:bioprj_PRJNA280774_Poecilia-latipinna/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA559876_Pomacanthus-maculosus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA379028_Porites-astreoides
rclone copy divdiv_datafiles:bioprj_PRJNA559876_Pomacanthus-maculosus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA379028_Porites-astreoides/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA544861_Porites-compressa
rclone copy divdiv_datafiles:bioprj_PRJNA379028_Porites-astreoides/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA544861_Porites-compressa/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA493660_Pygoscelis-adeliae
rclone copy divdiv_datafiles:bioprj_PRJNA544861_Porites-compressa/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA493660_Pygoscelis-adeliae/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA493660_Pygoscelis-antarcticus
rclone copy divdiv_datafiles:bioprj_PRJNA493660_Pygoscelis-adeliae/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA493660_Pygoscelis-antarcticus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA493660_Pygoscelis-papua
rclone copy divdiv_datafiles:bioprj_PRJNA493660_Pygoscelis-antarcticus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA493660_Pygoscelis-papua/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA573759_Pyroteuthis-margaritifera
rclone copy divdiv_datafiles:bioprj_PRJNA493660_Pygoscelis-papua/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA573759_Pyroteuthis-margaritifera/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA480308_Rhizophora-mangle
rclone copy divdiv_datafiles:bioprj_PRJNA573759_Pyroteuthis-margaritifera/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA480308_Rhizophora-mangle/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA553831_Robustosergia-robusta
rclone copy divdiv_datafiles:bioprj_PRJNA480308_Rhizophora-mangle/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA553831_Robustosergia-robusta/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA549138_Sargassum-muticum
rclone copy divdiv_datafiles:bioprj_PRJNA553831_Robustosergia-robusta/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA549138_Sargassum-muticum/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA253681_Scomber-scombrus
rclone copy divdiv_datafiles:bioprj_PRJNA549138_Sargassum-muticum/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA253681_Scomber-scombrus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA560239_Sebastes-diaconus
rclone copy divdiv_datafiles:bioprj_PRJNA253681_Scomber-scombrus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA560239_Sebastes-diaconus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA451040_Sebastes-paucispinis
rclone copy divdiv_datafiles:bioprj_PRJNA560239_Sebastes-diaconus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA451040_Sebastes-paucispinis/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA451040_Sebastes-pinniger
rclone copy divdiv_datafiles:bioprj_PRJNA451040_Sebastes-paucispinis/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA451040_Sebastes-pinniger/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA316872_Sebastes-rastrelliger
rclone copy divdiv_datafiles:bioprj_PRJNA451040_Sebastes-pinniger/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA316872_Sebastes-rastrelliger/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA451040_Sebastes-ruberrimus
rclone copy divdiv_datafiles:bioprj_PRJNA316872_Sebastes-rastrelliger/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA451040_Sebastes-ruberrimus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA359404_Sebastiscus-marmoratus
rclone copy divdiv_datafiles:bioprj_PRJNA451040_Sebastes-ruberrimus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA359404_Sebastiscus-marmoratus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA319656_Seriola-lalandi-dorsalis
rclone copy divdiv_datafiles:bioprj_PRJNA359404_Sebastiscus-marmoratus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA319656_Seriola-lalandi-dorsalis/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA385011_Siphamia-tubifer
rclone copy divdiv_datafiles:bioprj_PRJNA319656_Seriola-lalandi-dorsalis/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA385011_Siphamia-tubifer/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA286089_Sphyrna-tiburo
rclone copy divdiv_datafiles:bioprj_PRJNA385011_Siphamia-tubifer/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA286089_Sphyrna-tiburo/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA448430_Stegastes-beebei
rclone copy divdiv_datafiles:bioprj_PRJNA286089_Sphyrna-tiburo/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA448430_Stegastes-beebei/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA448430_Stegastes-leucorus
rclone copy divdiv_datafiles:bioprj_PRJNA448430_Stegastes-beebei/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA448430_Stegastes-leucorus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA361144_Stephanocoenia-intersepta
rclone copy divdiv_datafiles:bioprj_PRJNA448430_Stegastes-leucorus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA361144_Stephanocoenia-intersepta/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA655996_Strombus-pugilis
rclone copy divdiv_datafiles:bioprj_PRJNA361144_Stephanocoenia-intersepta/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA655996_Strombus-pugilis/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA354496_Symphodus-melops
rclone copy divdiv_datafiles:bioprj_PRJNA655996_Strombus-pugilis/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA354496_Symphodus-melops/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA361214_Symphodus-tinca
rclone copy divdiv_datafiles:bioprj_PRJNA354496_Symphodus-melops/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA361214_Symphodus-tinca/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA553831_Systellaspis-debilis
rclone copy divdiv_datafiles:bioprj_PRJNA361214_Symphodus-tinca/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA553831_Systellaspis-debilis/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA432036_Thunnus-thynnus
rclone copy divdiv_datafiles:bioprj_PRJNA553831_Systellaspis-debilis/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA432036_Thunnus-thynnus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA360511_Trachidermus-fasciatus
rclone copy divdiv_datafiles:bioprj_PRJNA432036_Thunnus-thynnus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA360511_Trachidermus-fasciatus/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA564451_Uria-aalge
rclone copy divdiv_datafiles:bioprj_PRJNA360511_Trachidermus-fasciatus/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA564451_Uria-aalge/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA564451_Uria-lomvia
rclone copy divdiv_datafiles:bioprj_PRJNA564451_Uria-aalge/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA564451_Uria-lomvia/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA573759_Vampyroteuthis-infernalis
rclone copy divdiv_datafiles:bioprj_PRJNA564451_Uria-lomvia/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA573759_Vampyroteuthis-infernalis/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA348572_Zalophus-wollebaeki
rclone copy divdiv_datafiles:bioprj_PRJNA573759_Vampyroteuthis-infernalis/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA348572_Zalophus-wollebaeki/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA449377_Zebrasoma-flavescens
rclone copy divdiv_datafiles:bioprj_PRJNA348572_Zalophus-wollebaeki/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA449377_Zebrasoma-flavescens/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJNA449377_Zebrasoma-scopas
rclone copy divdiv_datafiles:bioprj_PRJNA449377_Zebrasoma-flavescens/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJNA449377_Zebrasoma-scopas/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
echo starting dataset bioprj_PRJDB7819_Zenarchopterus-dunckeri
rclone copy divdiv_datafiles:bioprj_PRJNA449377_Zebrasoma-scopas/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"
rclone copy divdiv_datafiles:bioprj_PRJDB7819_Zenarchopterus-dunckeri/gendiv_data /Users/rachel/ALL_r80_gendiv_data --include "*WMfitwishart*"
rclone copy divdiv_datafiles:bioprj_PRJDB7819_Zenarchopterus-dunckeri/r80_outputs /Users/rachel/ALL_r80_popgen_data --include "*popgenstats.0.5*"



