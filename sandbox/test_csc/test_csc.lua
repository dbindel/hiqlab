A = CSCMatrix(4,4, 0,0);
B = CSCMatrix(4,4, 0,0);
P = CSCMatrix(4,2, 0,0);

A:load_triple({1, 2, 2, 3, 3, 4, 4},
              {1, 1, 2, 2, 3, 3, 4},
              {2, 1, 2, 1, 2, 1, 2},
              {1, 0, 0, 0, 0, 0, 0});
P:load_triple({3, 4}, 
              {1, 2}, 
              {1, 1})
B:load_triple({1, 1, 1, 1},
              {1, 2, 3, 4},
              {1, 2, 3, 4});

A:dump('A2.txt');
P:dump('P2.txt');

C = A:add(3, 2, B, nil)
C:dump('C.txt')

AT = A:transpose();
AT:dump('AT.txt');

ATA = AT:mul(A,nil);
ATA:dump('ATA.txt');

PT = P:transpose()
PTAP = A:mul(PT, P, nil)
PTAP:dump('PTAP.txt');

PTA = PT:mul(A,nil)
PTAP2 = PTA:mul(P,nil)
PTAP2:dump('PTAP2.txt');

C:delete()
AT:delete()
ATA:delete()
PT:delete()
PTA:delete()
PTAP:delete()
PTAP2:delete()
