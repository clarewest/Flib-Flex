# Convert the 8-category ss8 prediction to a 3-category prediction readable by Flib
# 3-category: H E C
# 8-category: H G I E B T S L
# H = H + G + I 
# E = E + B 
# C = T + S + L
awk '/^[^#]/ {printf "%4s %s %-3s %.3f %.3f %.3f\n", $1, $2, $3, $4+$5+$6, $7+$8, $9+$10+$11}' $1.ss8 > $1.fasta.ss 
