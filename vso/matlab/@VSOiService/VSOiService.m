function obj = VSOiService

obj.endpoint = 'http://vso.nso.edu/cgi-bin/VSO/PROD/vsoi_wsdl.cgi';
obj.wsdl = 'file:///C:/Users/mike/proj/solar/trunk/vso/matlab/VSOi_rpc_literal.wsdl';

obj = class(obj,'VSOiService');

