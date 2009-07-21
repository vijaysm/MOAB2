#!/bin/sh

default_fcdefs=iBase_FCDefs.h
case "x$1" in
  x|x-*)
    echo "Usage: $0 <PREFIX> <INPUT_FILE> <OUTPUT_FILE> [<FCDEFS>=${default_fcdefs}]"
    ;;
esac

PFX="$1"
INPUT="$2"
OUTPUT="$3"
if test "x" != "x$4"; then
  FCDEFS="$4"
else
  FCDEFS="$default_fcdefs"
fi

if test "x" = "x$SED"; then
  SED=`which sed`
fi

EXPR="s/^[[:space:]]*void[[:space:]][[:space:]]*${PFX}_\([_a-zA-Z0-9][_a-zA-Z0-9]*\)[[:space:]]*(.*\$/${PFX}_\1/p"

echo "#include \"$FCDEFS\"" > $OUTPUT
echo '#ifdef FC_FUNC_' >> $OUTPUT
echo >> $OUTPUT
for func in `$SED -n "$EXPR" $INPUT`; do
  lower=`echo $func | tr '[:upper:]' '[:lower:]'`
  upper=`echo $func | tr '[:lower:]' '[:upper:]'`
  echo "#define $func FC_FUNC_( $lower, $upper )" >> $OUTPUT
done
echo >> $OUTPUT
echo "#endif" >> $OUTPUT
