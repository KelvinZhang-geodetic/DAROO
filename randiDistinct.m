function v2 = randiDistinct(nLimit, len, v1)
%RANDIDISTINCT ������ v1 ��ȫ���ظ��������������
%
% v2 = randiDistinct(nMax, len, v1)
%   nMax  ������ȡֵ���ޣ�1:nMax��
%   len   ������������ v2 �ĳ���
%   v1    ����������һ���Ѵ�������
%   v2    ���������� v1 �޽������������
   if numel(nLimit)==1,nLimit=1:nLimit;end
    pool = setdiff(nLimit(1):nLimit(2), v1);     % ʣ�����ֵ
    if numel(pool) < len
        error('ʣ�����ֵ���㣬�޷�����ָ�����ȵĲ��ظ�����');
    end
    v2 = pool(randperm(numel(pool), len));  % �����len��
    v2=v2(:)';
end

