#/bin/sh
#"sort a big file alphabetical order for first column and numerical for second column"


##help
usage="$(basename "$0") INPUTFILE [-h] -- program to sort a big bed or vcf, do it ALPHABETICALLY and output to stdout. Works faster by creating subfiles."

while getopts ':hs:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    s) seed=$OPTARG
       ;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
done



###program

infile=$*

if  [ -d temptosort ]; then rm -r temptosort ; fi
mkdir temptosort 
awk -F "\t" '{print > "temptosort/"$1".tosort"}' $infile
for i in temptosort/*.tosort; do  echo $i ; sort -k2,2n $i > $i.s ; done
cat temptosort/*.s | grep -v ^temptosort >  /dev/stdout	
rm -r temptosort