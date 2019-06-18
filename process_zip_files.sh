for file in *.zip
do
	direc="${file%.*}"
	mkdir $direc
	mv $file $direc
	cd $direc
	unzip -j $file
	cd ..
done
