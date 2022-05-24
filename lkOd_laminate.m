
function [KE,dKE]=lkOd(angle);

T=angle;
z =[(45668133223994449*cos(T)^2)/3298534883328 - (4629803850984265*cos(T)*sin(T))/3298534883328 + 21253860423849683/3298534883328,    (45668133223994449*sin(2*T))/8796093022208 + (4629803850984265*cos(T)^2)/4398046511104 + 28429337708992111/8796093022208,  (4629803850984265*cos(T)*sin(T))/3298534883328 - (45668133223994449*cos(T)^2)/3298534883328 + 3540495746132599/3298534883328,    (45668133223994449*sin(2*T))/8796093022208 - (4629803850984265*sin(T)^2)/4398046511104 - 28429337708992111/8796093022208,   (4629803850984265*cos(T)*sin(T))/6597069766656 - (45668133223994449*cos(T)^2)/6597069766656 - 21253860423849683/6597069766656,    (4629803850984265*sin(T)^2)/4398046511104 - (45668133223994449*sin(2*T))/8796093022208 - 37688945410960641/8796093022208, (45668133223994449*cos(T)^2)/6597069766656 - (4629803850984265*cos(T)*sin(T))/6597069766656 - 28334851916114881/6597069766656,    (4629803850984265*sin(T)^2)/4398046511104 - (45668133223994449*sin(2*T))/8796093022208 + 28429337708992111/8796093022208;       
    (45668133223994449*sin(2*T))/8796093022208 + (4629803850984265*cos(T)^2)/4398046511104 + 28429337708992111/8796093022208, (4629803850984265*cos(T)*sin(T))/3298534883328 - (45668133223994449*cos(T)^2)/3298534883328 + 5576832803987011/274877906944, 37688945410960641/8796093022208 - (45668133223994449*cos(T)*sin(T))/4398046511104 - (4629803850984265*cos(T)^2)/4398046511104,  (4629803850984265*cos(T)*sin(T))/6597069766656 - (45668133223994449*cos(T)^2)/6597069766656 + 361110027247491/137438953472, - (4629803850984265*cos(T)^2)/4398046511104 - (45668133223994449*cos(T)*sin(T))/4398046511104 - 28429337708992111/8796093022208, (45668133223994449*cos(T)^2)/6597069766656 - (4629803850984265*cos(T)*sin(T))/6597069766656 - 5576832803987011/549755813888,      (45668133223994449*sin(2*T))/8796093022208 - (4629803850984265*sin(T)^2)/4398046511104 - 28429337708992111/8796093022208, (45668133223994449*cos(T)^2)/3298534883328 - (4629803850984265*cos(T)*sin(T))/3298534883328 - 7021272912976975/549755813888; 
    (4629803850984265*cos(T)*sin(T))/3298534883328 - (45668133223994449*cos(T)^2)/3298534883328 + 3540495746132599/3298534883328,    (4629803850984265*sin(T)^2)/4398046511104 - (45668133223994449*sin(2*T))/8796093022208 + 28429337708992111/8796093022208, (45668133223994449*cos(T)^2)/3298534883328 - (4629803850984265*cos(T)*sin(T))/3298534883328 + 21253860423849683/3298534883328,    (4629803850984265*sin(T)^2)/4398046511104 - (45668133223994449*sin(2*T))/8796093022208 - 37688945410960641/8796093022208,   (45668133223994449*cos(T)^2)/6597069766656 - (4629803850984265*cos(T)*sin(T))/6597069766656 - 28334851916114881/6597069766656,    (45668133223994449*sin(2*T))/8796093022208 - (4629803850984265*sin(T)^2)/4398046511104 - 28429337708992111/8796093022208, (4629803850984265*cos(T)*sin(T))/6597069766656 - (45668133223994449*cos(T)^2)/6597069766656 - 21253860423849683/6597069766656,    (45668133223994449*sin(2*T))/8796093022208 + (4629803850984265*cos(T)^2)/4398046511104 + 28429337708992111/8796093022208;
    (4629803850984265*cos(T)^2)/4398046511104 + (45668133223994449*cos(T)*sin(T))/4398046511104 - 37688945410960641/8796093022208,  (4629803850984265*cos(T)*sin(T))/6597069766656 - (45668133223994449*cos(T)^2)/6597069766656 + 361110027247491/137438953472,      (4629803850984265*sin(T)^2)/4398046511104 - (45668133223994449*sin(2*T))/8796093022208 - 37688945410960641/8796093022208, (4629803850984265*cos(T)*sin(T))/3298534883328 - (45668133223994449*cos(T)^2)/3298534883328 + 5576832803987011/274877906944,        (4629803850984265*sin(T)^2)/4398046511104 - (45668133223994449*sin(2*T))/8796093022208 + 28429337708992111/8796093022208, (45668133223994449*cos(T)^2)/3298534883328 - (4629803850984265*cos(T)*sin(T))/3298534883328 - 7021272912976975/549755813888, (4629803850984265*cos(T)^2)/4398046511104 + (45668133223994449*cos(T)*sin(T))/4398046511104 + 28429337708992111/8796093022208, (45668133223994449*cos(T)^2)/6597069766656 - (4629803850984265*cos(T)*sin(T))/6597069766656 - 5576832803987011/549755813888;
    (4629803850984265*cos(T)*sin(T))/6597069766656 - (45668133223994449*cos(T)^2)/6597069766656 - 21253860423849683/6597069766656,    (4629803850984265*sin(T)^2)/4398046511104 - (45668133223994449*sin(2*T))/8796093022208 - 37688945410960641/8796093022208, (45668133223994449*cos(T)^2)/6597069766656 - (4629803850984265*cos(T)*sin(T))/6597069766656 - 28334851916114881/6597069766656,    (4629803850984265*sin(T)^2)/4398046511104 - (45668133223994449*sin(2*T))/8796093022208 + 28429337708992111/8796093022208,   (45668133223994449*cos(T)^2)/3298534883328 - (4629803850984265*cos(T)*sin(T))/3298534883328 + 21253860423849683/3298534883328,    (45668133223994449*sin(2*T))/8796093022208 + (4629803850984265*cos(T)^2)/4398046511104 + 28429337708992111/8796093022208,  (4629803850984265*cos(T)*sin(T))/3298534883328 - (45668133223994449*cos(T)^2)/3298534883328 + 3540495746132599/3298534883328,    (45668133223994449*sin(2*T))/8796093022208 - (4629803850984265*sin(T)^2)/4398046511104 - 28429337708992111/8796093022208;
    - (4629803850984265*cos(T)^2)/4398046511104 - (45668133223994449*cos(T)*sin(T))/4398046511104 - 28429337708992111/8796093022208, (45668133223994449*cos(T)^2)/6597069766656 - (4629803850984265*cos(T)*sin(T))/6597069766656 - 5576832803987011/549755813888,      (45668133223994449*sin(2*T))/8796093022208 - (4629803850984265*sin(T)^2)/4398046511104 - 28429337708992111/8796093022208, (45668133223994449*cos(T)^2)/3298534883328 - (4629803850984265*cos(T)*sin(T))/3298534883328 - 7021272912976975/549755813888,        (45668133223994449*sin(2*T))/8796093022208 + (4629803850984265*cos(T)^2)/4398046511104 + 28429337708992111/8796093022208, (4629803850984265*cos(T)*sin(T))/3298534883328 - (45668133223994449*cos(T)^2)/3298534883328 + 5576832803987011/274877906944, 37688945410960641/8796093022208 - (45668133223994449*cos(T)*sin(T))/4398046511104 - (4629803850984265*cos(T)^2)/4398046511104,  (4629803850984265*cos(T)*sin(T))/6597069766656 - (45668133223994449*cos(T)^2)/6597069766656 + 361110027247491/137438953472;
    (45668133223994449*cos(T)^2)/6597069766656 - (4629803850984265*cos(T)*sin(T))/6597069766656 - 28334851916114881/6597069766656,    (45668133223994449*sin(2*T))/8796093022208 - (4629803850984265*sin(T)^2)/4398046511104 - 28429337708992111/8796093022208, (4629803850984265*cos(T)*sin(T))/6597069766656 - (45668133223994449*cos(T)^2)/6597069766656 - 21253860423849683/6597069766656,    (45668133223994449*sin(2*T))/8796093022208 + (4629803850984265*cos(T)^2)/4398046511104 + 28429337708992111/8796093022208,    (4629803850984265*cos(T)*sin(T))/3298534883328 - (45668133223994449*cos(T)^2)/3298534883328 + 3540495746132599/3298534883328,    (4629803850984265*sin(T)^2)/4398046511104 - (45668133223994449*sin(2*T))/8796093022208 + 28429337708992111/8796093022208, (45668133223994449*cos(T)^2)/3298534883328 - (4629803850984265*cos(T)*sin(T))/3298534883328 + 21253860423849683/3298534883328,    (4629803850984265*sin(T)^2)/4398046511104 - (45668133223994449*sin(2*T))/8796093022208 - 37688945410960641/8796093022208;
    (4629803850984265*sin(T)^2)/4398046511104 - (45668133223994449*sin(2*T))/8796093022208 + 28429337708992111/8796093022208, (45668133223994449*cos(T)^2)/3298534883328 - (4629803850984265*cos(T)*sin(T))/3298534883328 - 7021272912976975/549755813888, (4629803850984265*cos(T)^2)/4398046511104 + (45668133223994449*cos(T)*sin(T))/4398046511104 + 28429337708992111/8796093022208, (45668133223994449*cos(T)^2)/6597069766656 - (4629803850984265*cos(T)*sin(T))/6597069766656 - 5576832803987011/549755813888,   (4629803850984265*cos(T)^2)/4398046511104 + (45668133223994449*cos(T)*sin(T))/4398046511104 - 37688945410960641/8796093022208,  (4629803850984265*cos(T)*sin(T))/6597069766656 - (45668133223994449*cos(T)^2)/6597069766656 + 361110027247491/137438953472,      (4629803850984265*sin(T)^2)/4398046511104 - (45668133223994449*sin(2*T))/8796093022208 - 37688945410960641/8796093022208, (4629803850984265*cos(T)*sin(T))/3298534883328 - (45668133223994449*cos(T)^2)/3298534883328 + 5576832803987011/274877906944];
 
dw =[ (4629803850984265*sin(T)^2)/3298534883328 - (4629803850984265*cos(T)^2)/3298534883328 - (45668133223994449*cos(T)*sin(T))/1649267441664,                                             (45668133223994449*cos(2*T))/4398046511104 - (4629803850984265*cos(T)*sin(T))/2199023255552,  (4629803850984265*cos(T)^2)/3298534883328 - (4629803850984265*sin(T)^2)/3298534883328 + (45668133223994449*cos(T)*sin(T))/1649267441664,                                             (45668133223994449*cos(2*T))/4398046511104 - (4629803850984265*cos(T)*sin(T))/2199023255552,  (4629803850984265*cos(T)^2)/6597069766656 - (4629803850984265*sin(T)^2)/6597069766656 + (45668133223994449*cos(T)*sin(T))/3298534883328,                                             (4629803850984265*cos(T)*sin(T))/2199023255552 - (45668133223994449*cos(2*T))/4398046511104,  (4629803850984265*sin(T)^2)/6597069766656 - (4629803850984265*cos(T)^2)/6597069766656 - (45668133223994449*cos(T)*sin(T))/3298534883328,                                             (4629803850984265*cos(T)*sin(T))/2199023255552 - (45668133223994449*cos(2*T))/4398046511104;
    (45668133223994449*cos(2*T))/4398046511104 - (4629803850984265*cos(T)*sin(T))/2199023255552, (4629803850984265*cos(T)^2)/3298534883328 - (4629803850984265*sin(T)^2)/3298534883328 + (45668133223994449*cos(T)*sin(T))/1649267441664, (45668133223994449*sin(T)^2)/4398046511104 - (45668133223994449*cos(T)^2)/4398046511104 + (4629803850984265*cos(T)*sin(T))/2199023255552, (4629803850984265*cos(T)^2)/6597069766656 - (4629803850984265*sin(T)^2)/6597069766656 + (45668133223994449*cos(T)*sin(T))/3298534883328, (45668133223994449*sin(T)^2)/4398046511104 - (45668133223994449*cos(T)^2)/4398046511104 + (4629803850984265*cos(T)*sin(T))/2199023255552, (4629803850984265*sin(T)^2)/6597069766656 - (4629803850984265*cos(T)^2)/6597069766656 - (45668133223994449*cos(T)*sin(T))/3298534883328,                                              (45668133223994449*cos(2*T))/4398046511104 - (4629803850984265*cos(T)*sin(T))/2199023255552, (4629803850984265*sin(T)^2)/3298534883328 - (4629803850984265*cos(T)^2)/3298534883328 - (45668133223994449*cos(T)*sin(T))/1649267441664;
    (4629803850984265*cos(T)^2)/3298534883328 - (4629803850984265*sin(T)^2)/3298534883328 + (45668133223994449*cos(T)*sin(T))/1649267441664,                                             (4629803850984265*cos(T)*sin(T))/2199023255552 - (45668133223994449*cos(2*T))/4398046511104,  (4629803850984265*sin(T)^2)/3298534883328 - (4629803850984265*cos(T)^2)/3298534883328 - (45668133223994449*cos(T)*sin(T))/1649267441664,                                             (4629803850984265*cos(T)*sin(T))/2199023255552 - (45668133223994449*cos(2*T))/4398046511104,  (4629803850984265*sin(T)^2)/6597069766656 - (4629803850984265*cos(T)^2)/6597069766656 - (45668133223994449*cos(T)*sin(T))/3298534883328,                                             (45668133223994449*cos(2*T))/4398046511104 - (4629803850984265*cos(T)*sin(T))/2199023255552,  (4629803850984265*cos(T)^2)/6597069766656 - (4629803850984265*sin(T)^2)/6597069766656 + (45668133223994449*cos(T)*sin(T))/3298534883328,                                             (45668133223994449*cos(2*T))/4398046511104 - (4629803850984265*cos(T)*sin(T))/2199023255552;
    (45668133223994449*cos(T)^2)/4398046511104 - (45668133223994449*sin(T)^2)/4398046511104 - (4629803850984265*cos(T)*sin(T))/2199023255552, (4629803850984265*cos(T)^2)/6597069766656 - (4629803850984265*sin(T)^2)/6597069766656 + (45668133223994449*cos(T)*sin(T))/3298534883328,                                              (4629803850984265*cos(T)*sin(T))/2199023255552 - (45668133223994449*cos(2*T))/4398046511104, (4629803850984265*cos(T)^2)/3298534883328 - (4629803850984265*sin(T)^2)/3298534883328 + (45668133223994449*cos(T)*sin(T))/1649267441664,                                              (4629803850984265*cos(T)*sin(T))/2199023255552 - (45668133223994449*cos(2*T))/4398046511104, (4629803850984265*sin(T)^2)/3298534883328 - (4629803850984265*cos(T)^2)/3298534883328 - (45668133223994449*cos(T)*sin(T))/1649267441664, (45668133223994449*cos(T)^2)/4398046511104 - (45668133223994449*sin(T)^2)/4398046511104 - (4629803850984265*cos(T)*sin(T))/2199023255552, (4629803850984265*sin(T)^2)/6597069766656 - (4629803850984265*cos(T)^2)/6597069766656 - (45668133223994449*cos(T)*sin(T))/3298534883328;
    (4629803850984265*cos(T)^2)/6597069766656 - (4629803850984265*sin(T)^2)/6597069766656 + (45668133223994449*cos(T)*sin(T))/3298534883328,                                             (4629803850984265*cos(T)*sin(T))/2199023255552 - (45668133223994449*cos(2*T))/4398046511104,  (4629803850984265*sin(T)^2)/6597069766656 - (4629803850984265*cos(T)^2)/6597069766656 - (45668133223994449*cos(T)*sin(T))/3298534883328,                                             (4629803850984265*cos(T)*sin(T))/2199023255552 - (45668133223994449*cos(2*T))/4398046511104,  (4629803850984265*sin(T)^2)/3298534883328 - (4629803850984265*cos(T)^2)/3298534883328 - (45668133223994449*cos(T)*sin(T))/1649267441664,                                             (45668133223994449*cos(2*T))/4398046511104 - (4629803850984265*cos(T)*sin(T))/2199023255552,  (4629803850984265*cos(T)^2)/3298534883328 - (4629803850984265*sin(T)^2)/3298534883328 + (45668133223994449*cos(T)*sin(T))/1649267441664,                                             (45668133223994449*cos(2*T))/4398046511104 - (4629803850984265*cos(T)*sin(T))/2199023255552;
    (45668133223994449*sin(T)^2)/4398046511104 - (45668133223994449*cos(T)^2)/4398046511104 + (4629803850984265*cos(T)*sin(T))/2199023255552, (4629803850984265*sin(T)^2)/6597069766656 - (4629803850984265*cos(T)^2)/6597069766656 - (45668133223994449*cos(T)*sin(T))/3298534883328,                                              (45668133223994449*cos(2*T))/4398046511104 - (4629803850984265*cos(T)*sin(T))/2199023255552, (4629803850984265*sin(T)^2)/3298534883328 - (4629803850984265*cos(T)^2)/3298534883328 - (45668133223994449*cos(T)*sin(T))/1649267441664,                                              (45668133223994449*cos(2*T))/4398046511104 - (4629803850984265*cos(T)*sin(T))/2199023255552, (4629803850984265*cos(T)^2)/3298534883328 - (4629803850984265*sin(T)^2)/3298534883328 + (45668133223994449*cos(T)*sin(T))/1649267441664, (45668133223994449*sin(T)^2)/4398046511104 - (45668133223994449*cos(T)^2)/4398046511104 + (4629803850984265*cos(T)*sin(T))/2199023255552, (4629803850984265*cos(T)^2)/6597069766656 - (4629803850984265*sin(T)^2)/6597069766656 + (45668133223994449*cos(T)*sin(T))/3298534883328;
    (4629803850984265*sin(T)^2)/6597069766656 - (4629803850984265*cos(T)^2)/6597069766656 - (45668133223994449*cos(T)*sin(T))/3298534883328,                                             (45668133223994449*cos(2*T))/4398046511104 - (4629803850984265*cos(T)*sin(T))/2199023255552,  (4629803850984265*cos(T)^2)/6597069766656 - (4629803850984265*sin(T)^2)/6597069766656 + (45668133223994449*cos(T)*sin(T))/3298534883328,                                             (45668133223994449*cos(2*T))/4398046511104 - (4629803850984265*cos(T)*sin(T))/2199023255552,  (4629803850984265*cos(T)^2)/3298534883328 - (4629803850984265*sin(T)^2)/3298534883328 + (45668133223994449*cos(T)*sin(T))/1649267441664,                                             (4629803850984265*cos(T)*sin(T))/2199023255552 - (45668133223994449*cos(2*T))/4398046511104,  (4629803850984265*sin(T)^2)/3298534883328 - (4629803850984265*cos(T)^2)/3298534883328 - (45668133223994449*cos(T)*sin(T))/1649267441664,                                             (4629803850984265*cos(T)*sin(T))/2199023255552 - (45668133223994449*cos(2*T))/4398046511104
    (4629803850984265*cos(T)*sin(T))/2199023255552 - (45668133223994449*cos(2*T))/4398046511104, (4629803850984265*sin(T)^2)/3298534883328 - (4629803850984265*cos(T)^2)/3298534883328 - (45668133223994449*cos(T)*sin(T))/1649267441664, (45668133223994449*cos(T)^2)/4398046511104 - (45668133223994449*sin(T)^2)/4398046511104 - (4629803850984265*cos(T)*sin(T))/2199023255552, (4629803850984265*sin(T)^2)/6597069766656 - (4629803850984265*cos(T)^2)/6597069766656 - (45668133223994449*cos(T)*sin(T))/3298534883328, (45668133223994449*cos(T)^2)/4398046511104 - (45668133223994449*sin(T)^2)/4398046511104 - (4629803850984265*cos(T)*sin(T))/2199023255552, (4629803850984265*cos(T)^2)/6597069766656 - (4629803850984265*sin(T)^2)/6597069766656 + (45668133223994449*cos(T)*sin(T))/3298534883328,                                              (4629803850984265*cos(T)*sin(T))/2199023255552 - (45668133223994449*cos(2*T))/4398046511104, (4629803850984265*cos(T)^2)/3298534883328 - (4629803850984265*sin(T)^2)/3298534883328 + (45668133223994449*cos(T)*sin(T))/1649267441664];


w = double(z);
KE=w;
dKE=double(dw);
end