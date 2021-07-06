% findInterpFrames_verifyScript.m

if_1 = findInterpFrames([nan,0,1,0]);
if(all(if_1(:) == [0,0,1,0].'))
    disp('Test 1 Passed');
else
    disp('Test 1 Failed');
end

if_2 = findInterpFrames([nan,0,4,0,nan]);
if(all(if_2(:) == [0,0,1,0,0].'))
    disp('Test 2 Passed');
else
    disp('Test 2 Failed');
end

if_3 = findInterpFrames([4,0,1,1,1,0]);
if(all(if_3(:) == [0,0,1,1,1,0].'))
    disp('Test 3 Passed');
else
    disp('Test 3 Failed');
end

if_4 = findInterpFrames([nan,0,1,1,nan]);
if(all(if_4(:) == [0,0,0,0,0].'))
    disp('Test 4 Passed');
else
    disp('Test 4 Failed');
end

if_5 = findInterpFrames([nan,0,1,2,nan]);
if(all(if_5(:) == [0,0,0,0,0].'))
    disp('Test 5 Passed');
else
    disp('Test 5 Failed');
end

if_6 = findInterpFrames([nan,0,1,0,0,1,0]);
if(all(if_6(:) == [0,0,1,0,0,1,0].'))
    disp('Test 6 Passed');
else
    disp('Test 6 Failed');
end

if_7 = findInterpFrames([nan,0,1,1]);
if(all(if_7(:) == [0,0,0,0].'))
    disp('Test 7 Passed');
else
    disp('Test 7 Failed');
end




