#!/bin/bash
releaseToTest=$1

unpackDirectory='/tmp/Thesisvorlage/' # where can we unpack the files for testing?
texMasterfile='Thesis'

echo "Verifying Release Package ${releaseToTest}..."

if [ $# -eq 0 ]
then
    echo "!!! usage is ./verifyReleasePackage thePackageThatShouldBeChecked.zip"
    exit 1
fi


#check if file exists
if [ ! -f "$releaseToTest" ]; then
    echo "!!! file ${releaseToTest} not found."
    exit 2
fi
echo '[x] file exists'

#check if file is valid zip-archive
zip -T -q "${releaseToTest}"
zipProblem=$?
if [ ! $zipProblem -eq 0 ]
then
	echo "!!! ${releaseToTest}  is not a valid zip-archive"
	exit 3
fi
echo '[x] file is valid zip-archive'

#unpack archive
echo 'unpacking archive'
rm -rf "${unpackDirectory}"
mkdir "${unpackDirectory}" #-p: proceed without error if directory already exists, create parent directories if required
unzip -q -d "${unpackDirectory}" "${releaseToTest}"

#compile to pdf
cd "${unpackDirectory}" #step into unpacked directory so latex can find the included files

#initial run to generate the .aux file for bibtex
echo 'initial pdflatex run to generate .aux file for bibtex'
pdflatex -interaction=nonstopmode -output-directory="${unpackDirectory}" "${texMasterfile}.tex" > /dev/null
resultLatexInitial=$?

#bibtex for literature
echo 'bibtex'
bibtex "${texMasterfile}.aux" > /dev/null
resultBibtex=$?

#first proper run (errors not expected, but would be irrelevant)
echo 'first run'
pdflatex -interaction=nonstopmode -output-directory="${unpackDirectory}" "${texMasterfile}.tex" > /dev/null
resultLatexFirst=$?

#second proper run (for table of contents etc.)
echo 'second run'
pdflatex -interaction=nonstopmode -output-directory="${unpackDirectory}" "${texMasterfile}.tex" > /dev/null
resultLatexSecond=$?
echo 'compilation sequence finished'

#clean up
echo 'cleaning up temporary files'
rm -rf "${unpackDirectory}"

#show/return test result
if [ ! $resultLatexSecond -eq 0 ]
then
	echo '!!! latex file has errors'
	exit 4
fi
echo '[x] pdf can be compiled from LaTeX without errors'
exit 0
