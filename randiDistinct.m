function v2 = randiDistinct(nLimit, len, v1)
%RANDIDISTINCT 生成与 v1 完全无重复的随机整数数组
%
% v2 = randiDistinct(nMax, len, v1)
%   nMax  标量，取值上限（1:nMax）
%   len   标量，欲生成 v2 的长度
%   v1    行向量，第一组已存在数据
%   v2    行向量，与 v1 无交集的随机整数
   if numel(nLimit)==1,nLimit=1:nLimit;end
    pool = setdiff(nLimit(1):nLimit(2), v1);     % 剩余可用值
    if numel(pool) < len
        error('剩余可用值不足，无法生成指定长度的不重复数组');
    end
    v2 = pool(randperm(numel(pool), len));  % 随机抽len个
    v2=v2(:)';
end

