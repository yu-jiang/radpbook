para
P0=lyap((A-B*K)',Q+K'*R*K)
K1=inv(R)*B'*P0
P1=lyap((A-B*K1)',Q+K1'*R*K1)
K2=inv(R)*B'*P1
P2=lyap((A-B*K2)',Q+K2'*R*K2)
K3=inv(R)*B'*P2