#!/bin/bash

for dir in "$@"; do

    dir=`basename $dir`;

    dulist=(`du -s $dir`);

    #usedsectors=`perl -e "print ${dulist[0]} * 1.1 * 1024"`;
    usedsectors=`perl -e "print ${dulist[0]} * 2.1"`;

    hdiutil create -sectors $usedsectors -fs HFS+ -volname $dir.dmg $dir;

    hdiutil mount $dir.dmg;

    ditto -rsrcFork $dir /Volumes/$dir.dmg;

    hdiutil unmount /Volumes/$dir.dmg;

done;
