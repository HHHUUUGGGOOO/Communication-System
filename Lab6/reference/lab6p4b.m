len = 10^5;
N0 = 1;
BER = [];
for i = 0:2:20
    binary_data = floor(rand(1,len)*2);
    impulse_response = [1,0,1;1,1,1];
    encoded_data = conv_enc(binary_data, impulse_response);
    Eb = 10^(i/20)*N0;
    d = 2*Eb^(1/2);
    M = 2;
    sym_seq = symbol_mapper(encoded_data, M, d, 'PSK');
    for i = 1:length(sym_seq)
        ni = normrnd(0,N0);
        nq = normrnd(0,N0);
        sym_seq(i,1) = sym_seq(i,1)+ni;
        sym_seq(i,2) = sym_seq(i,2)+nq;
    end
    decoded_data = conv_dec2(sym_seq, impulse_response, d);
    err = sum(binary_data ~= decoded_data)/len
    BER = [BER err];
end
semilogy(0:2:20, BER)

%1
function proper_seq = generate_proper_seq(M)
    if M == 2
        proper_seq = [0;1];
    else
        seq = generate_proper_seq(M/2);
        seq_flip = flip(seq);
        x1 = zeros(M/2, 1);
        x2 = x1+1;
        seq1 = [x1 seq];
        seq2 = [x2 seq_flip];
        proper_seq = [seq1;seq2];
    end
end

function proper_seq2 = generate_proper_seq2(M)
    seq = generate_proper_seq(M^(1/2));
    proper_seq2 = [];
    for i = 1:M
        a = floor((i-1)/M^(1/2))+1;
        b = i-(a-1)*M^(1/2);
        proper_seq2 = [proper_seq2; [seq(a,:) seq(b,:)]];
    end
end

function sym_seq = symbol_mapper(bin_seq, M, d, name)
    l = fix(log2(M));
    if name == 'PAM'
        proper_seq = generate_proper_seq(M);
        sym_seq = [];
        idx = 1;
        while idx <= length(bin_seq)
            for m = 1:M
                if bin_seq(idx:idx+l-1) == proper_seq(m,:)
                    sym_seq = [sym_seq; [(2*m-1-M)*d/2 0]];
                end
            end
            idx = idx+l;
        end
    end
    if name == 'PSK'
        proper_seq = generate_proper_seq(M);
        sym_seq = [];
        idx = 1;
        while idx <= length(bin_seq)
            for m = 1:M
                if bin_seq(idx:idx+l-1) == proper_seq(m,:)
                    sym_seq = [sym_seq;[d/2/sin(pi/M)*cos(2*pi/M*(m-1)), d/2/sin(pi/M)*sin(2*pi/M*(m-1))]];
                end
            end
            idx = idx+l;
        end
    end
    if name == 'QAM'
        proper_seq2 = generate_proper_seq2(M);
        sym_seq = [];
        idx = 1;
        while idx <= length(bin_seq)
            for m = 1:M
                if bin_seq(idx:idx+l-1) == proper_seq2(m,:)
                    a = floor((m-1)/M^(1/2))+1;
                    b = m-(a-1)*M^(1/2);
                    sym_seq = [sym_seq;[(2*a-1-M^(1/2))*d/2 (2*b-1-M^(1/2))*d/2]];
                end
            end
            idx = idx+l;
        end
    end
end

%2
function bin_seq = MD_symbol_demapper(sym_seq, M, d, name)
    if name == 'QAM'
        proper_seq2 = generate_proper_seq2(M);
    else
        proper_seq2 = generate_proper_seq(M);
    end
    proper_seq3 = reshape(transpose(proper_seq2), [1 M*log2(M)]);
    sym_seq2 = symbol_mapper(proper_seq3, M, d, name);
    bin_seq = [];
    for i = 1:length(sym_seq)
        for j = 1:length(sym_seq2)
            if j == 1
                dmin = ((sym_seq(i,1)-sym_seq2(j,1))^2+(sym_seq(i,2)-sym_seq2(j,2))^2)^(1/2);
                c = j;
            else
                if dmin > ((sym_seq(i,1)-sym_seq2(j,1))^2+(sym_seq(i,2)-sym_seq2(j,2))^2)^(1/2)
                    dmin = ((sym_seq(i,1)-sym_seq2(j,1))^2+(sym_seq(i,2)-sym_seq2(j,2))^2)^(1/2);
                    c = j;
                end
            end
        end
        bin_seq = [bin_seq proper_seq2(c,:)];
    end
end

%functions from lab5
%2a
function encoded_data = conv_enc(binary_data, impulse_response)
    b = [0, 0, binary_data];
    encoded_data = [];
    for i = 1:length(binary_data)
        for j = 1:length(impulse_response(:,1))
            encoded_data = [encoded_data, mod(dot(flip(b(i:i+2)), impulse_response(j,:)), 2)];
        end
    end
end

%2b
function decoded_data = conv_dec2(binary_data, impulse_response, d)
    K = length(impulse_response(:,1));
    pathlen = length(binary_data)/K;
    statenum = 2^(length(impulse_response(1,:))-1);
    change = zeros(pathlen+1, statenum);
    V = zeros(pathlen+1, statenum);
    s = zeros(pathlen+1, statenum);
    s(1,1) = 1;
    all_coded = FSM(impulse_response);
    all_coded = -(all_coded-1/2)*d;
    for i = 2:length(V)
        for j = 1:statenum
            if s(i-1,j) ~= 0
                s1 = floor((j+1)/2);
                s2 = s1+statenum/2;
                p1 = s1*2-mod(j,2);
                p2 = s2*2-mod(j,2);
                dist1 = (sum((binary_data((i-2)*K+1:(i-1)*K) - all_coded(p1,:)).^(2)))^(1/2);
                if (V(i,s1) == 0 && change(i,s1) == 0)|| V(i,s1) > V(i-1,j)+dist1
                    V(i,s1) = V(i-1,j)+dist1;
                    s(i,s1) = j;
                    change(i,s1) = change(i,s1)+1;
                end
                dist2 = (sum((binary_data((i-2)*K+1:(i-1)*K) - all_coded(p2,:)).^(2)))^(1/2);
                if (V(i,s2) == 0 && change(i,s2) == 0) || V(i,s2) > V(i-1,j)+dist2
                    V(i,s2) = V(i-1,j)+dist2;
                    s(i,s2) = j;
                    change(i,s2) = change(i,s2)+1;
                end
            end
        end
    end
    s(1,1) = 0;
    state = 1;
    decoded_data = [0];
    for i = 1:pathlen-1
        if floor((state+1)/2) == 1
            decoded_data = [0, decoded_data];
        else
            decoded_data = [1, decoded_data];
        end
        state = s(pathlen+1-i, state);
    end
end
function all_coded = FSM(impulse_response)
    d = (0:2^length(impulse_response(1,:))-1)';
    b = de2bi(d);
    all_coded = zeros(length(b), length(impulse_response(:,1)));
    for i = 1:length(b)
        for j = 1:length(impulse_response(:,1))
            all_coded(i,j) = mod(dot(flip(b(i,:)), impulse_response(j,:)), 2);
        end
    end
end

%2c
function Y = mbsc(X, p)
    Y = [];
    p_array = rand(1,length(X));
    for i = 1:length(X)
        if p_array(i) < p
            Y(i) = mod(X(i)+1, 2);
        else
            Y(i) = X(i);
        end
    end
end