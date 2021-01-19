s_1 = {'s0', 's1', 's2', 's3', 's4'};
p_1 = [0.26, 0.25, 0.20, 0.15, 0.14];

% 2-(a) Dictionary
symbols = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'};
prob = [0.28, 0.21, 0.2, 0.12, 0.075, 0.06, 0.05, 0.005];
dict = huffman_dict( symbols, prob );
disp(dict)

% 2-(b) Encode
sym_seq = { 'g', 'a', 'c', 'a', 'b' };
bin_seq = huffman_enc(sym_seq, dict);
disp(bin_seq)

% 2-(c) Decode
sym_seq = huffman_dec(bin_seq, dict);
disp(sym_seq)

function d = huffman_dict(s, p)
    % Variables
    d = {};
    s_new = {};
    p_new = [];
    s_2 = s;
    p_2 = p;
    % Construct a huff-tree
    while (true)
        sum = 0;
        new_key = '';
        % Combine the 1st & 2nd min, replace them with their combination
        for i = 1:2
            ind = find(p_2==min(min(p_2)));
            new_key = [s_2{1, ind}(1, :), new_key];
            sum = sum + p_2(ind);
            s_2(ind) = [];
            p_2(ind) = [];
        end
        s_2 = [s_2, new_key];
        p_2 = [p_2, sum];
        s_new = [s_new, new_key];
        p_new = [p_new, sum];
        if (sum >= 1)
            break;
        end
    end
    % Create a tree
    temp_tree = {};
    tree = {};
    for i = 1:length(s_new)
        % Determine left / right nodes
        all_s = [s, s_new];
        all_p = [p, p_new];
       % Left / Right node default
       left = [];
       right = [];
       comp_p = p_new(i);
       for k = 1:length(all_s)
           alpha = all_p(k);
           for h = 1:length(all_s)
               beta = all_p(h);
               len = length(all_s{1, k}(1, :));
               if (alpha + beta == comp_p)
                   if (all_s{1, k}(1, :) == s_new{1, i}(1, 1:len))
                       left = [left, k];
                       right = [right, h];
                   end
               end
           end
       end
        temp_tree = {s_new{1, i}(1, :), left(1), right(1), i+length(s)};
        tree = [tree;temp_tree];
    end
    t = {};
    comb = {};
    all = {};
    for i = 1:length(s)
       left = [];
       right = [];
       t = {s{1, i}(1, :), left, right, i};
       comb = [comb;t];
    end
    all = [comb;tree];
    % Add the leaves
    for i = 1:length(s)
        % Left / Right node default
        left = [];
        right = [];
        enc = '';
        now = i; % index
        % Encoding
        while (true)
            for j = 1:length(all)
                if (all{j, 2} == now) % left node
                    enc = ['0', enc];
                    now = all{j, 4};
                    break;
                end
                if (all{j, 3} == now) % right node
                    enc = ['1', enc];
                    now = all{j, 4};
                    break;
                end
            end
            if (now == length(all))
                break;
            end
        end
        temp = {s{1, i}(1, :), p(i), left, right, enc};
        d = [d;temp];
    end
    % Add huff-tree other nodes
    for i = 1:length(s_new)
       left = [];
       right = [];
       enc = '';
       now = i+length(s); % index
       % Encoding
       while (true)
           for j = 1:length(all)
               if (all{j, 2} == now) % left node
                   enc = ['0', enc];
                   now = all{j, 4};
                   break;
               end
               if (all{j, 3} == now) % right node
                   enc = ['1', enc];
                   now = all{j, 4};
                   break;
               end
           end
           if (now == length(all))
               break;
           end
       end
       for j = 1:length(tree)
           if (length(s_new{1, i}(1, :)) == length(tree{j, 1}))
               if (s_new{1, i}(1, :) == tree{j, 1})
                   left = [left, tree{j, 2}];
                   right = [right, tree{j, 3}];
               end
           end
       end
       temp = {s_new{1, i}(1, :), p_new(i), left(1), right(1), enc};
       d = [d;temp];
    end
end

function enc = huffman_enc(ss, dic)
    enc = '';
    for i = 1:length(ss)
        for j = 1:length(dic)
            if (length(ss{1, i}) == length(dic{j, 1}))
                if (ss{1, i} == dic{j, 1})
                    enc = [enc, dic{j, 5}];
                end
            end
        end
    end
end

function dec = huffman_dec(b, dic)
    len = 0;
    dec = {};
    for i = 1:length(dic)
        if (length(dic{i, 1}) == 1)
            len = len+1;
        end
    end
    temp = length(b);
    while (temp > 0)
        for j = 1:len
            if (length(b) >= length(dic{j, 5}))
                if (dic{j, 5} == b(1:length(dic{j, 5})))
                    dec = [dec, dic{j, 1}];
                    b(1:length(dic{j, 5})) = [];
                    temp = length(b);
                    break;
                end
            end
        end
    end
end