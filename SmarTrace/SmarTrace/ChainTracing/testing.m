function testing()

test = struct
test = tester(test)
disp(test.one)
disp(test.two)
end

function test = tester(test)
test.one = magic(1)
test.two = magic(2)
end