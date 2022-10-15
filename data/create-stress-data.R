require(phytools)

parameters <- read.csv('wals/cldf/parameters.csv')
languages <- read.csv('wals/cldf/languages.csv')
values <- read.csv('wals/cldf/values.csv')
codes <- read.csv('wals/cldf/codes.csv')

merged <- merge(codes,values,by.x='ID',by.y='Code_ID')
stress1 <- merged[merged$Parameter_ID.x=='14A',c('Language_ID','Name')]
stress2 <- merged[merged$Parameter_ID.x=='15A',c('Language_ID','Name')]
stress <- droplevels(merge(stress1,stress2,by='Language_ID',all.x=T,all.y=T))

stress.bin <- cbind(to.matrix(stress$Name.x,seq=levels(stress$Name.x)),to.matrix(stress$Name.y,seq=levels(stress$Name.y)))
stress <- merge(stress,languages,by.x='Language_ID',by.y='ID')
#stress[stress$Language_ID=='jak',]$Glottocode <- 'popt1235'
stress.bin[stress.bin[,c("Combined: Right-edge and unbounded")]==1,c('Right-oriented: One of the last three','Unbounded: Stress can be anywhere')] <- 1
stress.bin <- stress.bin[,-which(colnames(stress.bin) %in% c('No fixed stress','Combined: Right-edge and unbounded','Fixed stress (no weight-sensitivity)'))]
stress.bin <- data.frame(stress.bin)
stress.bin$lang <- as.character(stress$Glottocode)
stress.bin$lang[stress.bin$lang==''] <- 'popt1235'

stress.bin[stress.bin$lang=='atam1239',]$lang <- 'urad1239'
stress.bin[stress.bin$lang=='awad1243',]$lang <- 'bhoj1244'
#stress.bin[stress.bin$lang=='bala1315',]$lang <- ''
stress.bin[stress.bin$lang=='band1339',]$lang <- 'midd1357'
stress.bin[stress.bin$lang=='barc1235',]$lang <- 'sout2826'
#stress.bin[stress.bin$lang=='basq1250',]$lang <- ''
#stress.bin[stress.bin$lang=='bawm1236',]$lang <- ''
#stress.bin[stress.bin$lang=='goro1270',]$lang <- ''
stress.bin[stress.bin$lang=='hala1252',]$lang <- 'jehh1245'
#stress.bin[stress.bin$lang=='hass1238',]$lang <- ''
#stress.bin[stress.bin$lang=='hija1235',]$lang <- ''
#stress.bin[stress.bin$lang=='juma1249',]$lang <- ''
stress.bin[stress.bin$lang=='kenu1236',]$lang <- 'kenu1243'
stress.bin[stress.bin$lang=='kewa1250',]$lang <- 'west2599'
stress.bin[stress.bin$lang=='kola1285',]$lang <- 'ujir1237'
stress.bin[stress.bin$lang=='latv1249',]$lang <- 'stan1325'
stress.bin[stress.bin$lang=='mayk1239',]$lang <- 'mayi1235'
#stress.bin[stress.bin$lang=='napu1241',]$lang <- ''
#stress.bin[stress.bin$lang=='narr1259',]$lang <- ''
stress.bin[stress.bin$lang=='nepa1252',]$lang <- 'nepa1254'
#stress.bin[stress.bin$lang=='ngur1261',]$lang <- ''
stress.bin[stress.bin$lang=='norw1259',]$lang <- 'norw1258'
stress.bin[stress.bin$lang=='nuau1240',]$lang <- 'sout2895'
#stress.bin[stress.bin$lang=='nyam1277',]$lang <- ''
stress.bin[stress.bin$lang=='osse1243',]$lang <- 'iron1242'
#stress.bin[stress.bin$lang=='paci1278',]$lang <- ''
stress.bin[stress.bin$lang=='plai1258',]$lang <- 'wood1236'
stress.bin[stress.bin$lang=='pola1255',]$lang <- 'kash1274'
#stress.bin[stress.bin$lang=='saka1289',]$lang <- ''
stress.bin[stress.bin$lang=='slav1253',]$lang <- 'sout2959'
stress.bin[stress.bin$lang=='sorb1249',]$lang <- 'lowe1385'
stress.bin[stress.bin$lang=='sout1528',]$lang <- 'serb1264'
stress.bin[stress.bin$lang=='stan1306',]$lang <- 'indo1316'
#stress.bin[stress.bin$lang=='ston1242',]$lang <- ''
#stress.bin[stress.bin$lang=='tamn1235',]$lang <- ''
#stress.bin[stress.bin$lang=='thay1248',]$lang <- ''
#stress.bin[stress.bin$lang=='tuka1247',]$lang <- ''
#stress.bin[stress.bin$lang=='wahg1249',]$lang <- ''
#stress.bin[stress.bin$lang=='weri1253',]$lang <- ''
#stress.bin[stress.bin$lang=='west2443',]$lang <- ''
#stress.bin[stress.bin$lang=='yine1238',]$lang <- ''
#stress.bin[stress.bin$lang=='yiry1247',]$lang <- ''

stress.bin <- aggregate(.~lang,stress.bin,FUN=max)
row.names(stress.bin) <- stress.bin$lang
stress.bin <- stress.bin[,2:ncol(stress.bin)]

colnames(stress.bin)[7:12] <- c('LeftEdge','LeftOriented','NotPredictable','RightEdge','RightOriented','Unbounded')

langs <- read.csv('asjp/languages.csv')

langs$treeID <- paste(langs$classification_wals,langs$Name,sep='.')

write.csv(file='stressMtx.csv',stress.bin,quote=F,row.names=rownames(stress.bin))
