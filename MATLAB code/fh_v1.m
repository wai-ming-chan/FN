clc; clear all;
rng(12345)

Ns = 3;
Nh = 8;
Nu = 9;
Nj = 1;

X_jam = randi(Ns, 1, Nh); % Given hopping pattern of Jammer
X_k = zeros(Nu,Nh);
FT_map = zeros(Ns,Nh);
FT_map2 = zeros(Ns,Nh);

%fprintf("total combination: %d\n", Ns^Nu);
O1_values=zeros(1,Ns^Nu);
O2_values=zeros(1,Ns^Nu);

% Case 1: minimizing hamming distance  
for t=1:Nh
  dist_max = -1; %Nu*(Nu-1);
  max_id = -1;

  for perm_id = 1:Ns^Nu 
    for idx_ue_k = 1:Nu
      xk_tmp(idx_ue_k) = mod( floor(perm_id/ Ns^(idx_ue_k-1)), Ns) + 1; 
    end

    dist_tmp = 0;
    for idx_ue1=1:Nu
      for idx_ue2=1:Nu
        if (idx_ue1 ~= idx_ue2 )
          dist_tmp = dist_tmp + dist_H( xk_tmp(idx_ue1), xk_tmp(idx_ue2) );
        end
      end

      dist_tmp = dist_tmp + dist_H(xk_tmp(idx_ue1), X_jam(t));
    end
    
    if (dist_tmp > dist_max)
      dist_max = dist_tmp;
      max_id = perm_id;
    end
    O1_values(perm_id) = dist_tmp;
  end

  % found optimal pattern 
  for idx_ue_k =1:Nu
    X_k(idx_ue_k, t) = mod( floor(max_id/ Ns^(idx_ue_k-1)), Ns) + 1; 
  end

  for k=1:Nu
    FT_map( X_k(k,t), t) = FT_map( X_k(k,t), t) + 1;
  end
  FT_map( X_jam(t), t) = FT_map( X_jam(t),t ) + 1;

end

% Case 2: maximizing success patterns
X_k2 = zeros(Nu, Nh);
for t=1:Nh
  succ_max = -1;
  max_id = -1;

  for perm_id = 1:Ns^Nu
    for idx_ue_k = 1:Nu
      xk_tmp(idx_ue_k) = mod( floor(perm_id/ Ns^(idx_ue_k-1)), Ns) + 1;
    end

    succ_tmp = 0;
    for f=1:Ns
      in1 = length(find(X_jam(t) == f));
      in2 = length(find( xk_tmp == f ));

      succ_tmp = succ_tmp + length(find( (1+in1)*in2 == 1 ));
    end
    
    if (succ_tmp > succ_max)
      succ_max = succ_tmp;
      max_id = perm_id;
    end
    O2_values(perm_id) = succ_tmp;
  end

  % found optimal pattern 
  for idx_ue_k = 1:Nu
    X_k2(idx_ue_k, t) = mod( floor(max_id/ Ns^(idx_ue_k-1)), Ns ) + 1;
  end

  for k=1:Nu
    FT_map2( X_k2(k,t), t) = FT_map2( X_k2(k,t), t) + 1;
  end
  FT_map2( X_jam(t), t) = FT_map2( X_jam(t),t) + 1;
end

fprintf("Jammer pattern:\n")
disp(X_jam)

fprintf("===== max hamming distance =====\n")
fprintf("UE patterns:\n")
disp(X_k)

fprintf("Freq-Time pattern:\n")
disp(FT_map)


fprintf("===== max success =====\n")
fprintf("UE patterns:\n")
disp(X_k2)
fprintf("Freq-Time pattern:\n")
disp(FT_map2)

%fprintf("===== O1-O2 =====\n")
%disp(O1_values-O2_values)

function [output] = dist_H(X1, X2)
  output = length( find( X1-X2 ~= 0 ) );
end
