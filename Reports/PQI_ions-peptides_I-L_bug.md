# PQI ions-peptides bug 
*where isoleucine is swaped with leucine*


*This bug is fixed with **cure_IL_NA** function*



###Case 1

####LSGGVAVLK

mascot:  
AT3G23990.1; AT2G33210.1

ions:  
AT3G23990.1; AT2G28000.1; AT2G33210.1  
AT2G33210.1; AT2G28000.1; AT3G23990.1

peptides:  
AT2G28000.1

####LSGGVAVIK
mascot:  
AT2G28000.1

ions:  
AT2G28000.1;AT2G33210.1;AT3G23990.1  
\#7516 \#17770 \#25076

Turn to LSGGVAVLK in peptides table.

----

###Case 2

####INLGVGAYR
ions:  
AT5G19550.1;AT4G31990.1;AT5G11520.1  

mascot:  
AT5G19550.1

peptides:  
AT4G31990.1

####LNLGVGAYR
mascot:  
AT4G31990.1

ions:  
AT4G31990.1;AT5G11520.1;AT5G19550.1  
\#9699 \#19654 \#32613 \#55918

Turn to INLGVGAYR in peptides table
