function MECH413Project
%----------------------------------------------------------
% Derviş Barış Ercüment
% Efe Şen
%----------------------------------------------------------

% Data definition
%---------------------------------------------------------- 
Coords=table2array(importdata("Coordinates.mat"));              %Importing coordinates.
Conn=table2array(importdata("Connectivity.mat"));               %Importing connectivity matrix.
E = 200000;                                                     %Young's modulus of material in [MPa].
v = 0.3;                                                        %Poisson's ratio of material.
t = 10;                                                         %Thickness of elements in [mm].
%----------------------------------------------------------

% Material Property Matrix
%----------------------------------------------------------
C = E/((1+v)*(1-2*v)) * [1-v v 0; v 1-v 0; 0 0 (1-2*v)/2];      %Material property matrix for plane strain.
%C= E/(1-v^2) * [1 v 0 ; v 1 0 ; 0 0 (1-v)/2];                  %Material property matrix for plane stress.
%----------------------------------------------------------

% Area calculation
%----------------------------------------------------------
for i=1:1:size(Conn,1)                                          %Loops through connectivity matrix to get nodes corresponding to each element.
        node1=Conn(i,1);                                        %Node numbers are obtained here.
        node2=Conn(i,2);
        node3=Conn(i,3);
        x1=Coords(node1,1);                                     %Coordinates of nodes are obtained here.
        x2=Coords(node2,1);
        x3=Coords(node3,1);
        y1=Coords(node1,2);
        y2=Coords(node2,2);
        y3=Coords(node3,2);
        Area(i,1)=0.5*det([1 x1 y1; 1 x2 y2; 1 x3 y3]);         %Areas of the elements are obtained here.
end
%----------------------------------------------------------

% Local stiffness matrices
%----------------------------------------------------------
for i=1:1:size(Conn,1)                                          %Loops through connectivity matrix to get nodes corresponding to each element.
        node1=Conn(i,1);                                        %Node numbers are obtained here.
        node2=Conn(i,2);
        node3=Conn(i,3);
        x1=Coords(node1,1);                                     %Coordinates of nodes are obtained here.
        x2=Coords(node2,1);
        x3=Coords(node3,1);
        y1=Coords(node1,2);
        y2=Coords(node2,2);
        y3=Coords(node3,2);
        b1=y2-y3;                                               %Entries of element-strain displacement matrix are calculated here.
        b2=y3-y1;
        b3=y1-y2;
        c1=x3-x2;
        c2=x1-x3;
        c3=x2-x1;
        B = 1/(2*Area(i,1)) * [b1 0 b2 0 b3 0; 0 c1 0 c2 0 c3; c1 b1 c2 b2 c3 b3];
        K{i}=transpose(B)*C*B*Area(i,1)*t;                      %Cell array of elemental stiffness matrices.
end

        KL=[K{1}];                                              %Define K local.

for i=2:1:size(Conn,1)
        KL=[KL;K{i}];                                           %All elemental stiffness matrices concatenated vertically.
end
%----------------------------------------------------------

% Assembly of global stiffness matrix
%----------------------------------------------------------
ElNum = size(Conn,1);
LNodeNum = size(Conn,2);
GNodeNum = 223;

KG=zeros(2*GNodeNum,2*GNodeNum);                                %Preallocates the memory for the global stiffness matrix.

selector=0;
for RowNumConn=1:1:ElNum                                        %Traverses rows of connectivity matrix.
    for ColNumConn=1:1:LNodeNum                                 %Traverses columns of connectivity matrix.
        indexi=Conn(RowNumConn,ColNumConn);                     %Obtains connectivity matrix entry to determine targeted rows of KL.
        for dummy=1:1:LNodeNum                                  %Traverses columns of connectivity matrix for each iteration of indexi.
            indexj=Conn(RowNumConn,dummy);                      %Obtains connectivity matrix entry to determine targeted columns of KL.
            KLRowNumX=(2*ColNumConn-1)+(2*LNodeNum*selector);   %Selects row number of entry of KL to be inserted into KG for degree of freedom x.
            KLRowNumY=(2*ColNumConn)+(2*LNodeNum*selector);     %Selects row number of entry of KL to be inserted into KG for degree of freedom y.
            KLColNumX=2*dummy-1;                                %Selects column number of entry of KL to be inserted into KG for degree of freedom x.
            KLColNumY=2*dummy;                                  %Selects column number of entry of KL to be inserted into KG for degree of freedom y.
            z=(indexi-1)*2+1;                                   %Row number of KG corresponding to selected row number of KL for degree of freedom y.
            q=(indexj-1)*2+1;                                   %Column number of KG corresponding to selected column number of KL for degree of freedom y.
            KG(z,q)=KG(z,q)+KL(KLRowNumX,KLColNumX);            %Positions of KG at which values of KL must be placed.
            KG(z,q+1)=KG(z,q+1)+KL(KLRowNumX,KLColNumY);
            KG(z+1,q)=KG(z+1,q)+KL(KLRowNumY,KLColNumX);
            KG(z+1,q+1)=KG(z+1,q+1)+KL(KLRowNumY,KLColNumY);
        end
    end
    selector=selector+1;
end
%----------------------------------------------------------

%Boundary Conditions
%----------------------------------------------------------
% Nodes 1-17 have fixed displacement.
% Node 73 has Fx = 2800 N and Fy = 0 N.
% Node 78 has Fx = 0 N and Fy = 5000 N.
%----------------------------------------------------------

%Condensed Equations for primary variables.
%----------------------------------------------------------
ZeroDisplacement = 1:1:34;                                      %Rows and columns that correspond to the fixed nodes.
CKG=KG;                                                         %Condensed K Global.
CKG(ZeroDisplacement,:)=[];                                     %Removal of rows corresponding to zero displacement.
CKG(:,ZeroDisplacement)=[];                                     %Removal of columns corresponding to zero displacement.

F=zeros(412,1);                                                 %Column vector of known forces.
F(111)=2800;
F(122)=5000;
%----------------------------------------------------------

%Primary Variables
%----------------------------------------------------------
U=inv(CKG)*F;                                                   %Solution of primary variables.

for i=1:2:411
    Ux(i,1)=U(i,1);                                             %Horizontal displacements.
end

for i=2:2:412
    Uy(i,1)=U(i,1);                                             %Vertical displacements.
end
MaxVerticalDisp=max(Uy);
MaxHorizontalDisp=max(Ux);

U1=zeros(34,1);
U=[U1;U];                                                       %Complete set of displacements.
%----------------------------------------------------------

%Element stresses
%----------------------------------------------------------
for i=1:1:ElNum
    for j=1:1:LNodeNum
        value=Conn(i,j);
        switch j
            case 1
                U2(1,1)=U(2*value-1,1);                         
                U2(2,1)=U(2*value,1);
            case 2
                U2(3,1)=U(2*value-1,1);
                U2(4,1)=U(2*value,1);
            case 3
                U2(5,1)=U(2*value-1,1);
                U2(6,1)=U(2*value,1);
            otherwise
        end
    end
    U3{i}=U2;
end

for i=1:1:size(Conn,1)                                          %Loops through connectivity matrix to get nodes corresponding to each element.
        node1=Conn(i,1);                                        %Node numbers are obtained here.
        node2=Conn(i,2);
        node3=Conn(i,3);
        x1=Coords(node1,1);                                     %Coordinates of nodes are obtained here.
        x2=Coords(node2,1);
        x3=Coords(node3,1);
        y1=Coords(node1,2);
        y2=Coords(node2,2);
        y3=Coords(node3,2);
        b1=y2-y3;                                               %Entries of element-strain displacement matrix are calculated here.
        b2=y3-y1;
        b3=y1-y2;
        c1=x3-x2;
        c2=x1-x3;
        c3=x2-x1;
        B = 1/(2*Area(i,1)) * [b1 0 b2 0 b3 0; 0 c1 0 c2 0 c3; c1 b1 c2 b2 c3 b3];
        Stresses{i}=C*B*U3{i};
end

for i=1:1:ElNum
    ElStress=Stresses{i};
    Strx=ElStress(1,1);
    Stry=ElStress(2,1);
    Strxy=ElStress(3,1);
    VonMises(i,1)=sqrt(((Strx)^2)-(Strx*Stry)+((Stry)^2)+3*(Strxy)^2);
end
[MaxStress,MaxStressElementNumber] = max(VonMises);

fprintf('Maximum x displacement is %d [mm]\n',MaxHorizontalDisp)
fprintf('Maximum y displacement is %d [mm]\n',MaxVerticalDisp)
fprintf('Maximum von Mises stress is %d [MPa]\n',MaxStress)
fprintf('Maximum von Mises stress occurs in element %d\n\n',MaxStressElementNumber)
%----------------------------------------------------------

end