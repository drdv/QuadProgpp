%-------------------------------------------------
% [x, info] = quadprogpp (H)

TEST_ID = '1';
N = 20;
[x, info] = quadprogpp(eye(N));
xref = zeros(N,1);

if ((abs(sum(x-xref)) < 1e-08) && (info.obj == 0) && (info.status == 0))
    printf(['Test ', TEST_ID, ' OK\n']);
else
    printf(['Test ', TEST_ID, ' FAIL\n']);
end
%-------------------------------------------------


%-------------------------------------------------
% [x, info] = quadprogpp (H, g)

TEST_ID = '2';
N = 20;
[x, info] = quadprogpp(eye(N), ones(N, 1));
xref = -1*ones(N,1);

if ((abs(sum(x-xref)) < 1e-08) && (info.obj == 0) && (info.status == 0))
    printf(['Test ', TEST_ID, ' OK\n']);
else
    printf(['Test ', TEST_ID, ' FAIL\n']);
end
%-------------------------------------------------


%-------------------------------------------------
% [x, info] = quadprogpp (H, g, A, b)

TEST_ID = '3';
N = 20;
[x, info] = quadprogpp(eye(N), ones(N, 1), [eye(4), zeros(4, N-4)], [1; 2; 3; 4]);
xref = [1   2   3   4  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1]';

if ((abs(sum(x-xref)) < 1e-08) && (info.obj == 0) && (info.status == 0))
    printf(['Test ', TEST_ID, ' OK\n']);
else
    printf(['Test ', TEST_ID, ' FAIL\n']);
end
%-------------------------------------------------


%-------------------------------------------------
% [x, info] = quadprogpp (H, g, A, b, lb, ub) {infeasible}

TEST_ID = '4';
N = 20;
[x, info] = quadprogpp( eye(N),
                        ones(N, 1),
                        [eye(4), zeros(4, N-4)],
                        [1; 2; 3; 4],
                        -100*ones(N,1),
                        zeros(N,1));

if (info.status > 0)
    printf(['Test ', TEST_ID, ' OK\n']);
else
    printf(['Test ', TEST_ID, ' FAIL\n']);
end
%-------------------------------------------------


%-------------------------------------------------
% [x, info] = quadprogpp (H, g, A, b, lb, ub)

TEST_ID = '5';
N = 20;
[x, info] = quadprogpp( eye(N),
                        ones(N, 1),
                        [eye(4), zeros(4, N-4)],
                        [1; 2; 3; 4],
                        [-100*ones(4,1); -0.5*ones(N-4, 1)],
                        [100*ones(4,1); 0.5*ones(N-4, 1)]);
xref = [1.0   2.0   3.0   4.0  -0.5  -0.5 -0.5  -0.5  -0.5  -0.5  -0.5  -0.5 -0.5  -0.5  -0.5  -0.5  -0.5  -0.5 -0.5  -0.5]';


if ((abs(sum(x-xref)) < 1e-08) && (info.obj == 0) && (info.status == 0))
    printf(['Test ', TEST_ID, ' OK\n']);
else
    printf(['Test ', TEST_ID, ' FAIL\n']);
end
%-------------------------------------------------


%-------------------------------------------------
% [x, info] = quadprogpp (H, g, A, b, lb, ub, Ain, lbin, ubin)

TEST_ID = '6';
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

if ((abs(sum(x-xref)) < 1e-08) && (info.obj == 0) && (info.status == 0))
    printf(['Test ', TEST_ID, ' OK\n']);
else
    printf(['Test ', TEST_ID, ' FAIL\n']);
end
%-------------------------------------------------
