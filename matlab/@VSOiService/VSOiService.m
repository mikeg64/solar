function obj = VSOiService

obj.endpoint = 'http://vso.nascom.nasa.gov/cgi-bin/VSOi_strict';
obj.wsdl = 'file:///C:/Users/mikeg/proj/solar/matlab/sdacVSOi_strict.wsdl';

obj = class(obj,'VSOiService');

