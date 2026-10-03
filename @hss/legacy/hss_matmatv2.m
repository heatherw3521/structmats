function H = hss_matmatv2(A,B)
%HSS_MATMAT mult of two hss matrices
if A.size(2) ~= B.size(1)
    error("Mismatched sizes for matrix multiplication.")
elseif A.levelcount ~= B.levelcount
    error("HSS matrices must have the same structure for efficient matmat multiplication.")
end
if A.isleaf
    H = A;
    H.D = A.D*B.D;
    H.size = size(H.D);
    return
end
g = bottom_up(A, B);
f = 0;
H = A;
H = up_bottom(A, B, g, f, H);
end

function g = bottom_up(A,B)
g = struct();
if A.A11.isleaf==0
    g.gl = bottom_up(A.A11, B.A11);
    g.gr = bottom_up(A.A22, B.A22);
    if A.isroot == 0
        lp = A.A12.Y * g.gl.YZ * B.A21.Z;
        rp = A.A21.Y * g.gr.YZ * B.A12.Z;
        g.YZ = lp + rp;
    end
else  % at parent of leaves
    % CORRECTED: gl corresponds to left child (A11, B11) interaction
    % The left child's off-diagonal blocks are A12 (upper-right) and B21 (lower-left)
    g.gl.YZ = A.A12.lrcomponent * A.A12.Y * B.A21.Z * B.A21.lrcomponent;
    
    % CORRECTED: gr corresponds to right child (A22, B22) interaction  
    % The right child's off-diagonal blocks are A21 (lower-left) and B12 (upper-right)
    g.gr.YZ = A.A21.lrcomponent * A.A21.Y * B.A12.Z * B.A12.lrcomponent;
    
    lp = A.A12.Y * g.gl.YZ * B.A21.Z;
    rp = A.A21.Y * g.gr.YZ * B.A12.Z;
    g.YZ = lp + rp;
end
end

function H = up_bottom(A,B,g,f,H)
if A.A11.isleaf == 0
    if A.isroot
        fl = g.gr.YZ;
        fr = g.gl.YZ;
    else
        fl = g.gr.YZ + A.A21.Z * f * B.A12.Y;
        fr = g.gl.YZ + A.A12.Z * f * B.A21.Y;
    end
    
    if A.isroot
        % At root: C12 = A11*B12 + A12*B22
        % Keep separate rank contributions from A and B
        H.A12.lrcomponent = blkdiag(A.A12.lrcomponent, B.A12.lrcomponent);
        H.A21.lrcomponent = blkdiag(A.A21.lrcomponent, B.A21.lrcomponent);
        
        % Update bases to match the block-diagonal lrcomponent structure
        % Z: [A's contribution | B's contribution] - horizontal concat
        H.A12.Z = [A.A12.Z, B.A12.Z];  % [m × k_A, m × k_B] → m × (k_A + k_B)
        H.A21.Z = [A.A21.Z, B.A21.Z];
        
        % Y: [A's contribution; B's contribution] - vertical concat (Y is transposed)
        H.A12.Y = [A.A12.Y; B.A12.Y];  % [(k_A × n); (k_B × n)] → (k_A + k_B) × n
        H.A21.Y = [A.A21.Y; B.A21.Y];
    else
        H.A12.lrcomponent = [A.A12.lrcomponent A.A12.Z * f * B.A12.Y; 
                             zeros(size(B.A12.lrcomponent,1),size(A.A12.lrcomponent,2)) B.A12.lrcomponent];
        H.A21.lrcomponent = [A.A21.lrcomponent A.A21.Z * f * B.A21.Y; 
                             zeros(size(B.A21.lrcomponent,1),size(A.A21.lrcomponent,2)) B.A21.lrcomponent];
        
        % Reference: C.Wl = [A.Wl zeros; B.B21'*g.gr.VU'*A.Wr B.Wl]
        % Wl' = A12.Y (transposed), so we build the untransposed version then transpose
        temp_Y12 = [A.A12.Y' zeros(size(A.A12.Y',1),size(B.A12.Y',2)); 
                    B.A21.lrcomponent' * g.gr.YZ' * A.A21.Y' B.A12.Y'];
        H.A12.Y = temp_Y12';
        
        temp_Y21 = [A.A21.Y' zeros(size(A.A21.Y',1),size(B.A21.Y',2)); 
                    B.A12.lrcomponent' * g.gl.YZ' * A.A12.Y' B.A21.Y'];
        H.A21.Y = temp_Y21';
        
        % Z updates (horizontal block structure)
        H.A21.Z = [A.A21.Z A.A21.lrcomponent*g.gr.YZ*B.A12.Z; 
                   zeros(size(B.A21.Z,1),size(A.A21.Z,2)) B.A21.Z];
        H.A12.Z = [A.A12.Z A.A12.lrcomponent*g.gl.YZ*B.A21.Z; 
                   zeros(size(B.A12.Z,1),size(A.A12.Z,2)) B.A12.Z];
    end
    
    H.A11 = up_bottom(A.A11, B.A11, g.gl, fl, H.A11);
    H.A22 = up_bottom(A.A22, B.A22, g.gr, fr, H.A22);
    
else  % at leaf level
    if A.isroot
        fl = g.gr.YZ;
        fr = g.gl.YZ;
    else
        fl = g.gr.YZ + A.A21.Z * f * B.A12.Y;
        fr = g.gl.YZ + A.A12.Z * f * B.A21.Y;
    end
    
    % CORRECTED: expand fl and fr through lrcomponent before adding to D
    H.A11.D = A.A11.D * B.A11.D + A.A12.Z * A.A12.lrcomponent * fl * B.A21.lrcomponent * B.A21.Y;
    H.A22.D = A.A22.D * B.A22.D + A.A21.Z * A.A21.lrcomponent * fr * B.A12.lrcomponent * B.A12.Y;
    
    % CORRECTED: proper basis updates
    % C12 = A11*B12 + A12*B22
    %     = A11*B12 + A12.Z * A12.lrc * A12.Y * B22
    % New basis: Z from A12, Y from B12, plus A11*B12 contribution
    H.A12.Z = [A.A12.Z, A.A11.D * B.A12.Z];  % horizontal: [10×k_A, 10×k_B] → 10×(k_A+k_B)
    H.A12.Y = [A.A12.Y * B.A22.D; B.A12.Y];  % VERTICAL: Y is transposed, [(k_A×10)', (k_B×10)'] → (k_A+k_B)×10
    H.A12.lrcomponent = [A.A12.lrcomponent, zeros(size(A.A12.lrcomponent,1), size(B.A12.lrcomponent,2)); 
                         zeros(size(B.A12.lrcomponent,1), size(A.A12.lrcomponent,2)), B.A12.lrcomponent];
    
    % C21 = A21*B11 + A22*B21
    H.A21.Z = [A.A21.Z, A.A22.D * B.A21.Z];  % horizontal: [10×k_A, 10×k_B] → 10×(k_A+k_B)
    H.A21.Y = [A.A21.Y * B.A11.D; B.A21.Y];  % VERTICAL: Y is transposed
    H.A21.lrcomponent = [A.A21.lrcomponent, zeros(size(A.A21.lrcomponent,1), size(B.A21.lrcomponent,2));
                         zeros(size(B.A21.lrcomponent,1), size(A.A21.lrcomponent,2)), B.A21.lrcomponent];
end
end