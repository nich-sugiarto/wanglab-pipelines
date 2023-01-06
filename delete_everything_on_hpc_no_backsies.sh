#!/bin/bash

echo "Why would you run this"
echo "This script will self destruct in"

for ((n=5;n>-1;n--)); do
	echo ${n}...
	sleep 1
done

echo "Deleting now..."
sleep 5
echo done

rm -- $0