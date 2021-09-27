mkdir w_21
mkdir w_25
mkdir w_41
mkdir reindeer

needle="PATH TO Needle exectuable."

$needle minimiser -k 21 -w 21 -o w_21/  "PATH_TO_DIR"/SRR*
$needle minimiser -k 21 -w 25 -o w_25/  "PATH_TO_DIR"/SRR*
$needle minimiser -k 21 -w 41 -o w_41/  "PATH_TO_DIR"/SRR*
