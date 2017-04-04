function dtree( a,b,c,d )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
clc;
train=importdata(a);
test=importdata(b);
option=c;
prun=str2num(d);
[m,n]=size(train);
[tm,tn]=size(test);
un_class=unique(train(:,n));
number_class=length(un_class);
examples=train;
default=distribution(examples);
attributes=zeros(n-1,1);
for i=1:n-1
    attributes(i,1)=i;
end
if strcmp(option,'optimized') || strcmp(option,'randomized')
node=1;
tree=dtl(1,node,examples,attributes,default,prun,option);

tot=0.00;
for i=1:tm
    classified=false;
root=tree;
    accuracy=0.00;
    while classified==false
        if root.best_attribute==-1
            classified=true;
            actual=test(i,tn);
            high=max(root.distribution(:,2));
            in=find(high==root.distribution(:,2));
            ind=0;
            competing=zeros(length(in),1);
            if length(in)>1
                
                for j=1:length(competing)
                    competing(j,1)=root.distribution(in(j,1),1);
                end
                ind=in(1,1);
            else
                ind=in(1,1);
            end
            predicted=root.distribution(ind,1);
            if predicted==actual
                accuracy=1.00;
            elseif predicted~=actual
                if ismember(predicted,competing)
                    accuracy=1/length(in);
                else
                accuracy=0.00;
                end
            end
        else
        attr=root.best_attribute;
        thres=root.best_threshold;
        if test(i,attr)<thres
            root=root.left_child;
        else
            root=root.right_child;
        end
        end
    end
    fprintf('ID=%5d, predicted=%3d, true=%3d, accuracy=%4.2f\n', i-1, predicted, actual, accuracy);
    tot=tot+accuracy;
    
end
tot=tot/tm;
fprintf('classification accuracy=%6.4f\n', tot);
end

if strcmp(option,'forest3')
node=1;
tree1=dtl(1,node,examples,attributes,default,prun,option);
node=1;
tree2=dtl(2,node,examples,attributes,default,prun,option);
node=1;
tree3=dtl(3,node,examples,attributes,default,prun,option);
end

if strcmp(option,'forest15')
node=1;
tree1=dtl(1,node,examples,attributes,default,prun,option);
node=1;
tree2=dtl(2,node,examples,attributes,default,prun,option);
node=1;
tree3=dtl(3,node,examples,attributes,default,prun,option);
node=1;
tree4=dtl(4,node,examples,attributes,default,prun,option);
node=1;
tree5=dtl(5,node,examples,attributes,default,prun,option);
node=1;
tree6=dtl(6,node,examples,attributes,default,prun,option);
node=1;
tree7=dtl(7,node,examples,attributes,default,prun,option);
node=1;
tree8=dtl(8,node,examples,attributes,default,prun,option);
node=1;
tree9=dtl(9,node,examples,attributes,default,prun,option);
node=1;
tree10=dtl(10,node,examples,attributes,default,prun,option);
node=1;
tree11=dtl(11,node,examples,attributes,default,prun,option);
node=1;
tree12=dtl(12,node,examples,attributes,default,prun,option);
node=1;
tree13=dtl(13,node,examples,attributes,default,prun,option);
node=1;
tree14=dtl(14,node,examples,attributes,default,prun,option);
node=1;
tree15=dtl(15,node,examples,attributes,default,prun,option);
end
end


function[tree]= dtl(tree_n,node,examples, attributes, default,pruning,opt)
[p,q]=size(unique(examples(:,end)));
[em,en]=size(examples);
option=opt;
if length(examples)<pruning
    tree=struct('node',node,'best_attribute',-1,'best_threshold',-1,'distribution',default);
    fprintf('tree=%2d,node=%3d, feature=%2d, thr=%6.2f, gain=%f\n',tree_n-1, node, -1, -1, -1);
elseif p==1 && q==1
    tree=struct('node',node,'best_attribute',-1,'best_threshold',-1,'distribution',distribution(examples));
    fprintf('tree=%2d,node=%3d, feature=%2d, thr=%6.2f, gain=%f\n',tree_n-1, node, -1, -1, -1);
else
    if strcmp(opt,'optimized')
     [best_attribute, best_threshold,gain,left_sum,right_sum] = choose_attribute(examples, attributes);
    elseif strcmp(opt,'randomized') || strcmp(opt,'forest15') || strcmp(opt,'forest3')
        [best_attribute, best_threshold,gain,left_sum,right_sum] = choose_attribute_ran(examples, attributes);
    end
     tree = struct('node',node,'best_attribute',best_attribute,'best_threshold',best_threshold,'distribution',distribution(examples));
     fprintf('tree=%2d,node=%3d, feature=%2d, thr=%6.2f, gain=%f\n',tree_n-1, node, best_attribute-1, best_threshold, gain);
     left_examples=zeros(left_sum,en);
     right_examples=zeros(right_sum,en);
     left_count=1;
     right_count=1;
     for i=1:em
         if examples(i,best_attribute)<best_threshold           
                 left_examples(left_count,:)=examples(i,:);
                 left_count=left_count+1;
         end
         if examples(i,best_attribute)>=best_threshold
                 right_examples(right_count,:)=examples(i,:);
                 right_count=right_count+1;
         end    
     end
     pruning_t=pruning;
     examples_left = left_examples;
     examples_right = right_examples;
     if length(examples_left)<pruning_t
         tree.left_child = dtl(tree_n,2*node,examples_left, attributes, default,pruning_t,option);     
     else
         tree.left_child = dtl(tree_n,2*node,examples_left, attributes, distribution(examples_left),pruning_t,option);
     end
     if length(examples_right)<pruning_t
         tree.right_child = dtl(tree_n,(2*node)+1,examples_right, attributes, default,pruning_t,option);     
     else
     tree.right_child = dtl(tree_n,(2*node)+1,examples_right, attributes, distribution(examples_right),pruning_t,option);
     end
end
end

function [best_attribute,best_threshold,gain,left_sum,right_sum] = choose_attribute(examples, attributes) 
 max_gain = -1;
 best_attribute = -1;
 best_threshold = -1;
 [natt,~]=size(attributes);
 for i=1:natt
     attribute_values = examples(:,i);
     L = min(attribute_values);
     M = max(attribute_values); 
     for k=1:50
         threshold = L + k*(M-L)/51;
         [gain,left_s,right_s] = info_gain(examples, attributes(i,1), threshold);
         if gain > max_gain
             max_gain = gain;
             best_attribute = attributes(i,1);
             best_threshold = threshold;
             left_sum=left_s;
             right_sum=right_s;
         end
     end
 end
end

function [best_attribute,best_threshold,gain,left_sum,right_sum] = choose_attribute_ran(examples, attributes) 
 max_gain = -1;
 best_attribute = -1;
 best_threshold = -1;
 natt=datasample(attributes,1);
     attribute_values = examples(:,natt);
     L = min(attribute_values);
     M = max(attribute_values); 
     for k=1:50
         threshold = L + k*(M-L)/51;
         [gain,left_s,right_s] = info_gain(examples, natt, threshold);
         if gain > max_gain
             max_gain = gain;
             best_attribute = natt;
             best_threshold = threshold;
             left_sum=left_s;
             right_sum=right_s;
         end
     end
end

function [gain,left_sum,right_sum]= info_gain(examples,A,threshold)
temp=examples(:,A);
[num,~]=size(unique(examples(:,end)));
un_class=unique(examples(:,end));
left=zeros(num,2);
right=zeros(num,2);
left(:,1)=un_class(:,1);
right(:,1)=un_class(:,1);
[s,~]=size(temp);
[em,en]=size(examples);
for i=1:s
    if temp (i,1)<threshold
        temp_class=examples(i,en);
        for j=1:num
            if left(j,1)==temp_class
                left(j,2)=left(j,2)+1;
            end
        end
    end
    if temp(i,1)>=threshold
        temp_class=examples(i,en);
        for j=1:num
            if right(j,1)==temp_class
                right(j,2)=right(j,2)+1;
            end
        end
    end    
end
root=zeros(num,2);
root(:,1)=un_class(:,1);
for i=1:num
    root(i,2)=left(i,2)+right(i,2);
end
root_sum=sum(root(:,2));
left_sum=sum(left(:,2));
right_sum=sum(right(:,2));
root_ent=0.0;
for i=1:num
    if root(i,2)~=0
        root_ent=root_ent-((root(i,2)/root_sum)*log2(root(i,2)/root_sum));
    end
end
left_ent=0.0;
for i=1:num
    if left(i,2)~=0
        left_ent=left_ent-((left(i,2)/left_sum)*log2(left(i,2)/left_sum));
    end
end
right_ent=0.0;
for i=1:num
    if right(i,2)~=0
        right_ent=right_ent-((right(i,2)/right_sum)*log2(right(i,2)/right_sum));
    end
end
gain=root_ent - ((left_sum/root_sum)*left_ent + (right_sum/root_sum)*right_ent );
end

function [dist]=distribution(examples)
[dm,dn]=size(examples);
classes=unique(examples(:,dn));
n=length(classes);
dist=zeros(n,2);
dist(:,1)=classes(:,1);
for i=1:dm
    temp_class=examples(i,dn);
    ind=find(temp_class==dist(:,1));
    dist(ind,2)=dist(ind,2)+1;
end
for i=1:n
    dist(i,2)=dist(i,2)/dm;
end
end

