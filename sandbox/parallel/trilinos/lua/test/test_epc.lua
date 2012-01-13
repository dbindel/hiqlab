-- Test Epetra_Vector_Complex class
numid = 5
numvec= 1

-- Testing Copy constructor
A = Epetra_MultiVector_Complex:create(numid,numvec)
A:Random()
normA = {}
B = Epetra_MultiVector_Complex:new(A)
normB = {}
A:Scale(0.0,0.0)
A:Norm2(normA)
B:Norm2(normB)
print('Test 1: Copy constructor')
for i = 1, numvec do
    print('NormA:',normA[i],'  NormB:',normB[i])
end
print(' ')
B:delete()

-- Testing View constructor
A:Random()
normA = {}
--B = Epetra_MultiVector_Complex:new(A)
B = Epetra_MultiVector_Complex:new(View, A, 0, 2)
normB1 = {}
normB2 = {}
B:Norm2(normB1)
A:Scale(0.0,0.0)
A:Norm2(normA)
B:Norm2(normB2)
print('Test 2: View constructor')
for i = 1, numvec do
    print('NormA:',normA[i],'  NormB1:',normB1[i], '  NormB2:',normB2[i])
end
print(' ')
B:delete()

-- Testing Scale1 function
A:Random()
normA1 = {}
A:Norm2(normA1)
A:Scale(1,0)
normA2 = {}
A:Norm2(normA2)
A:Scale(0,1)
normA3 = {}
A:Norm2(normA3)
print('Test 3: Scale function')
for i = 1, numvec do
    print('NormA1:',normA1[i],'  NormA2:',normA2[i], '  NormA3:',normA3[i])
end
print(' ')

-- Testing Scale2 function
A:Random()
B = Epetra_MultiVector_Complex:new(A)
B:Random()
B:Norm2(normB)
normA1 = {}
A:Norm2(normA1)
A:Scale(1,0,B)
normA2 = {}
A:Norm2(normA2)
A:Scale(0,1,B)
normA3 = {}
A:Norm2(normA3)
A:Scale(0,0,B)
normA4 = {}
A:Norm2(normA4)
print('Test 4: Add Scale Vector function')
for i = 1, numvec do
    print('NormB :',normB[i] , '  NormA1:',normA1[i], '  NormA2:',normA2[i])
    print('NormB :',normB[i] , '  NormA3:',normA3[i], '  NormA4:',normA4[i])
end
print(' ')

-- Testing Update1 function
A:Random()
B:Random()
A:Norm2(normA1)
A:Update(1, 0, A, 0, 0)
A:Norm2(normA2)
A:Update(0, 1, A, 0, 0)
A:Norm2(normA3)
A:Update(0, 0, A, 1, 0)
A:Norm2(normA4)
A:Update(0, 0, A, 0, 1)
normA5 = {}
A:Norm2(normA5)
A:Update(0, 0, A, 0, 0)
normA6 = {}
A:Norm2(normA6)
print('Test 5: Test Update1 function')
for i = 1, numvec do
    print('NormA1:',normA1[i], '  NormA2:',normA2[i])
    print('NormA3:',normA3[i], '  NormA4:',normA4[i])
    print('NormA5:',normA5[i], '  NormA6:',normA6[i])
end
print(' ')

-- Testing Update2 function
A:Random()
B:Random()
A:Norm2(normA1)
A:Update(1, 0, A, 0, 0, A, 0, 0)
A:Norm2(normA2)
A:Update(0, 1, A, 0, 0, A, 0, 0)
A:Norm2(normA3)
A:Update(0, 0, A, 1, 0, A, 0, 0)
A:Norm2(normA4)
A:Update(0, 0, A, 0, 1, A, 0, 0)
normA5 = {}
A:Norm2(normA5)
A:Update(0, 0, A, 0, 0, A, 1, 0)
normA6 = {}
A:Norm2(normA6)
A:Update(0, 0, A, 0, 0, A, 0, 1)
normA7 = {}
A:Norm2(normA7)
A:Update(0, 0, A, 0, 0, A, 0, 0)
normA8 = {}
A:Norm2(normA8)
print('Test 6: Test Update2 function')
for i = 1, numvec do
    print('NormA1:',normA1[i], '  NormA2:',normA2[i])
    print('NormA3:',normA3[i], '  NormA4:',normA4[i])
    print('NormA5:',normA5[i], '  NormA6:',normA6[i])
    print('NormA7:',normA7[i], '  NormA8:',normA8[i])
end
print(' ')
