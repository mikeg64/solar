loadatmos;

loadcutoffdata;

plot(height(1420:2048),pres(1420:2048),height(1420:2048),pfit1);
plot(height(1420:2048)./1e6,cutoff_chromos(height(1420:2048)),hc./1e6,cutoff,height(1324:1419)./1e6,cutoff_transition(height(1324:1419)));
plot(height(1420:2048)./1e6,cutoff_chromos(height(1420:2048)),hc./1e6,cutoff,height(1324:1419)./1e6,cutoff_transition(height(1324:1419)),height(1:1325)./1e6,cutoff_corona(height(1:1325)));
