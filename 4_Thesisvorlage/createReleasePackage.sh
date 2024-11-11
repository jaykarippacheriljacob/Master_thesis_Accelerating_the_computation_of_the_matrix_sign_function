#!/bin/bash
version=$(date +"%Y-%m-%d")
workingDirectory=$PWD
outfile="${workingDirectory}/Thesisvorlage_${version}.zip"
tempDirectoryVersionIntegration='/tmp/Thesisvorlage/CopyForAutomaticVersionIntegration/'
texMasterfile='Thesis'

filesNeeded="*.tex siunitx.cfg ./Medien/*.* ./Kapitel/*.* ./Verzeichnisse/*.* ./Vorlage/*.* Thesisbeispiel-FAQ-Tipps.pdf readme.txt"
#######################################

filesNeeded_filtered=""
filesNeeded_pathsInTempdirectory=""
declare -A filesNeededBySuffix # associative array
declare -A suffixesFound

for file in ${filesNeeded}
do
	fileSuffix=${file##*.}
	fileName=${file%.*}
	#printf "%s : %s\n" ${fileSuffix} ${fileName};

	case $fileSuffix in
		tex|bib|cfg|png|pdf|txt)
			#echo "$file accepted"
			if [[ -v suffixesFound[$fileSuffix] ]]; then
				((suffixesFound[$fileSuffix]++))
			else
				suffixesFound[$fileSuffix]=1
			fi
			
			filesNeededBySuffix["${fileSuffix}-${suffixesFound[$fileSuffix]}"]=$file

			filesNeeded_filtered+=" ${file}"
			filesNeeded_pathsInTempdirectory+=" ${tempDirectoryVersionIntegration}${file}"
			;;
		
		*)
#			echo "$file ignored"
			;;
	esac
done

#echo "Suffixes found:"
#for suffix in "${!suffixesFound[@]}"
#do
#	printf " - %s (%d)\n" $suffix ${suffixesFound[$suffix]}
#	for (( index = 1; index <= ${suffixesFound[$suffix]}; index++ ));	do
#		printf "      %s\n" ${filesNeededBySuffix["$suffix-$index"]}
#	done
#done

#echo $filesNeeded_pathsInTempdirectory
#echo "${filesNeeded_filtered}"
#######################################

#clean up temp directory
echo 'Clean up temporary files (in case any are left from a previous run)'
rm -rf "${tempDirectoryVersionIntegration}"
mkdir -p "${tempDirectoryVersionIntegration}"

echo 'Create copy of required files'
cp --parents ${filesNeeded_filtered} "${tempDirectoryVersionIntegration}"

echo 'Go to copied directory'
cd ${tempDirectoryVersionIntegration}

echo 'Add info header to files where possible'

declare -a infoHeader
infoHeader+=("Version ${version}")
while read -r line
do
	infoHeader+=("${line}")
done < "${workingDirectory}/releaseInformationToPrepend"

for suffix in "${!suffixesFound[@]}"
do
	case $suffix in
		tex | bib | cfg | txt)
			case $suffix in
				tex | bib | cfg)
					blockSymbol="%"
					;;
				txt)
					blockSymbol="*"
					;;
				*)
					printf "ERROR: Unkown Case should never occur here as all suffixes have been checked before\n"
					exit 1
			esac

			for (( index = 1; index <= ${suffixesFound[$suffix]}; index++ ));	do
				currentFile=${filesNeededBySuffix["$suffix-$index"]}
				#echo "   Processing $currentFile"
				fileContent=$(<$currentFile)
				{
					for line in "${infoHeader[@]}"
					do
						echo "${blockSymbol}${blockSymbol} ${line}"
					done
					for i in {1..80}
					do
						printf "%s" "${blockSymbol}"
					done
					printf "\n"
					echo ""
					echo "${fileContent}"
				}>$currentFile
			done
			;;

		*)
			##echo  "   i dont know how to add text to ${suffix}-files"
			;;
	esac
done

printf 'Generate example pdf and check if original files compile without errors\n'
printf '\tInit-Compile for .aux-Files\n'
pdflatex -interaction=nonstopmode "${texMasterfile}.tex" > /dev/null

printf '\tBibTeX\n'
bibtex Thesis.aux > /dev/null

printf '\tWe need two compilation cycles for correct example pdf\n'
printf '\t\tfirst run...'
pdflatex -interaction=nonstopmode Thesis.tex > /dev/null
printf 'done\n'
printf '\t\tsecond run...'
pdflatex -interaction=nonstopmode Thesis.tex > /dev/null
printf 'done\n'
printf '\trenaming example pdf for end user\n'
cp Thesis.pdf 'Thesisbeispiel-FAQ-Tipps.pdf'

printf "Creating zip archive..."
rm -f "${outfile}"
zip -q -9 "${outfile}" ${filesNeeded}
#  zip -> zip archiv erstellen
#  -9  -> Stärkste Komprimierung verwenden (0 bis 9)
# outfile -> Name der Zip-Datei
#  Rest-> Dateien, die einbezogen werden sollen
printf "done.\n\tCurrent release has been packed as ${outfile}\n"
printf "\n"
cd ${workingDirectory}
# Noch mal entpacken und kompilieren, um die Vollständigkeit zu verifizieren?"
./verifyReleasePackage.sh "${outfile}"
verificationResult=$?
if [ $verificationResult -eq 0 ]
then
	printf 'Package is VALID\n'
else
	printf 'Package is INVALID\n'
	printf 'Please check if the archive contains all necessary files.\n'
	printf 'If that is the case, check if there are errors in the LaTeX-code.\n'
	exit 1
fi
exit 0


#TODO : Release-Version in die .tex-Dateien integrieren (Kommentar am Anfang)
# dafür müsste man:
# - eine Kopie aller Dateien erzeugen
# - in allen Dateien am Anfang den Versionsblock hinzufügen
# - dann mit der Kopie den Release erstellen
# - die Kopie löschen

#TODO : Automatic numbering *a, *b, *c for same-day releases
