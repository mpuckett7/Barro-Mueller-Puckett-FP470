for nodes in 1 2 4 8; do 
    for tpn in 1 2 4; do
        grep -oP "Post collide time:\s*\K[0-9.+e-]+" comtime-$nodes-$tpn.out > raw-$nodes-$tpn.txt
    done
done