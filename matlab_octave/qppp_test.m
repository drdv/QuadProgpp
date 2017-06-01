TEST_ID = 1;

%-------------------------------------------------
% [x, info] = quadprogpp (H)

N = 20;
[x, info] = quadprogpp(eye(N));
xref = zeros(N,1);

if ((abs(sum(x-xref)) < 1e-08) && (info.status == 0))
    printf('Test %d OK\n', TEST_ID);
else
    printf('Test %d FAIL\n', TEST_ID);
end
TEST_ID = TEST_ID + 1;
%-------------------------------------------------


%-------------------------------------------------
% [x, info] = quadprogpp (H, g)

N = 20;
[x, info] = quadprogpp(eye(N), ones(N, 1));
xref = -1*ones(N,1);

if ((abs(sum(x-xref)) < 1e-08) && (info.status == 0))
    printf('Test %d OK\n', TEST_ID);
else
    printf('Test %d FAIL\n', TEST_ID);
end
TEST_ID = TEST_ID + 1;
%-------------------------------------------------


%-------------------------------------------------
% [x, info] = quadprogpp (H, g, A, b)

N = 20;
[x, info] = quadprogpp(eye(N), ones(N, 1), [eye(4), zeros(4, N-4)], [1; 2; 3; 4]);
xref = [1   2   3   4  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1]';

if ((abs(sum(x-xref)) < 1e-08) && (info.status == 0))
    printf('Test %d OK\n', TEST_ID);
else
    printf('Test %d FAIL\n', TEST_ID);
end
TEST_ID = TEST_ID + 1;
%-------------------------------------------------


%-------------------------------------------------
% [x, info] = quadprogpp (H, g, A, b, lb, ub) {infeasible}

N = 20;
[x, info] = quadprogpp( eye(N),
                        ones(N, 1),
                        [eye(4), zeros(4, N-4)],
                        [1; 2; 3; 4],
                        -100*ones(N,1),
                        zeros(N,1));

if (info.status > 0)
    printf('Test %d OK\n', TEST_ID);
else
    printf('Test %d FAIL\n', TEST_ID);
end
TEST_ID = TEST_ID + 1;
%-------------------------------------------------


%-------------------------------------------------
% [x, info] = quadprogpp (H, g, A, b, lb, ub)

N = 20;
[x, info] = quadprogpp( eye(N),
                        ones(N, 1),
                        [eye(4), zeros(4, N-4)],
                        [1; 2; 3; 4],
                        [-100*ones(4,1); -0.5*ones(N-4, 1)],
                        [100*ones(4,1); 0.5*ones(N-4, 1)]);
xref = [1.0   2.0   3.0   4.0  -0.5  -0.5 -0.5  -0.5  -0.5  -0.5  -0.5  -0.5 -0.5  -0.5  -0.5  -0.5  -0.5  -0.5 -0.5  -0.5]';


if ((abs(sum(x-xref)) < 1e-08) && (info.status == 0))
    printf('Test %d OK\n', TEST_ID);
else
    printf('Test %d FAIL\n', TEST_ID);
end
TEST_ID = TEST_ID + 1;
%-------------------------------------------------


%-------------------------------------------------
% [x, info] = quadprogpp (H, g, A, b, lb, ub, Ain, lbin, ubin)

N = 20;
[x, info] = quadprogpp( eye(N),
                        ones(N, 1),
                        [eye(4), zeros(4, N-4)],
                        [1; 2; 3; 4],
                        [-100*ones(4,1); -5*ones(N-4, 1)],
                        [100*ones(4,1); 0.5*ones(N-4, 1)],
                        ones(1,N),
                        [-1.5],
                        [1.5]);
xref = [1.0   2.0   3.0   4.0  -0.71875  -0.71875 -0.71875  -0.71875  -0.71875  -0.71875  -0.71875  -0.71875 -0.71875  -0.71875  -0.71875  -0.71875  -0.71875  -0.71875 -0.71875  -0.71875]';

if ((abs(sum(x-xref)) < 1e-08) && (info.status == 0))
    printf('Test %d OK\n', TEST_ID);
else
    printf('Test %d FAIL\n', TEST_ID);
end
TEST_ID = TEST_ID + 1;
%-------------------------------------------------


%-------------------------------------------------
% [x, info] = quadprogpp (H, g, A, b) {linear dependence}

N = 20;
[x, info] = quadprogpp(eye(N), ones(N, 1), [[eye(4), zeros(4, N-4)]; [1, zeros(1, N-1)]], [1; 2; 3; 4; 1]);

if (info.status > 0)
    printf('Test %d OK\n', TEST_ID);
else
    printf('Test %d FAIL\n', TEST_ID);
end
TEST_ID = TEST_ID + 1;
%-------------------------------------------------
