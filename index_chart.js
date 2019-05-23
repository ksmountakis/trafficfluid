
function chart_show(vm) {
	let conf = {
		bindto: "#chart",
		data: {},
		step: 1
	};
	conf.step = Math.max(2, Math.round(vm.sim.K / 400));
	if (vm.chart.general.type === "mdsd")
		chart_show_mdsd(vm, conf);
	if (vm.chart.general.type === "mmux")
		chart_show_mmux(vm, conf);
	if (vm.chart.general.type === "flow")
		chart_show_flow(vm, conf);
	if (vm.chart.general.type === "reg")
		chart_show_reg(vm, conf);
	if (vm.chart.general.type === "frame")
		chart_show_frame(vm, conf);
}

function chart_show_flow(vm, conf) {
	console.log("showing flow");
	
	var flow = ["flow"];
	for (var t = 0; t < vm.sim.K; t += conf.step)
		flow.push(vm.sim.flow[t]);
	
	conf.data.columns = [flow];
	c3.generate(conf);
}

function chart_show_frame(vm, conf) {
	console.log("showing frame");
	
	var degree_max = ["degree_max"];
	var degree_avg = ["degree_avg"];
	for (var t = 0; t < vm.sim.K; t += 1) {
		degree_max.push(vm.sim.frame[t].max);
		degree_avg.push(vm.sim.frame[t].avg);
	}
	
	conf.data.columns = [degree_max, degree_avg];
	c3.generate(conf);
}

function chart_show_reg(vm, conf) {
	console.log("showing reg");
	
	var reg = ["reg"];
	for (var t = 0; t < vm.sim.K; t += conf.step)
		reg.push(vm.sim.reg[t]);
	
	conf.data.columns = [reg];
	c3.generate(conf);
}

function chart_show_mmux(vm, conf) {
	console.log("showing mmux");
	
	var maxux = ["maximum_ux"];
	var minux = ["minimum_ux"];
	var mean_up = ["mean_above_ux"];
	var mean_dn = ["mean_below_ux"];

	for (var t = 0; t < vm.sim.K; t += 1) {

		var max_t = -1000;
		var min_t = +1000;
		var mean_up_t = 0;
		var mean_dn_t = 0;
		for (var i = 0; i < vm.sim.n; i++) {
			let ux = vm.sim.ux[i][t];
			max_t = Math.max(max_t, ux);
			min_t = Math.min(min_t, ux);
			mean_up_t += Math.max(0, ux);
			mean_dn_t += Math.min(0, ux);
		 }
		 
		 mean_up_t /= vm.sim.n;
		 mean_dn_t /= vm.sim.n;

		maxux.push(max_t);
		minux.push(min_t);
		mean_up.push(mean_up_t);
		mean_dn.push(mean_dn_t);
	}
	conf.data.columns = [maxux, minux, mean_up, mean_dn];
	c3.generate(conf);
}

function chart_show_mdsd(vm, conf) {
	console.log("showing mdsd");

	var mdsd_up_max = ["mdsd_above_max"];
	var mdsd_up = ["mdsd_above"];
	var mdsd_dn = ["mdsd_below"];
	var mdsd_bo = ["mdsd_both"];

	for (var t = 0; t < vm.sim.K; t += 1) {
		var mdsd_up_max_t = 0;
		var mdsd_up_t = 0;
		var mdsd_bo_t = 0;
		var mdsd_dn_t = 0;
		for (var i = 0; i < vm.sim.n; i++) {
			let dv = Math.round(10000*(vm.sim.vx[i][t] - vm.sim.vd[i])/(0.1 + vm.sim.vd[i]))/100;
			mdsd_up_max_t = Math.max(dv, mdsd_up_max_t);
			mdsd_up_t += Math.max(dv, 0);
			mdsd_dn_t += Math.min(dv, 0);
			mdsd_bo_t += dv;
		}
		mdsd_up_t /= vm.sim.n;
		mdsd_dn_t /= vm.sim.n;
		mdsd_bo_t /= vm.sim.n;

		mdsd_up_max.push(mdsd_up_max_t);
		mdsd_up.push(mdsd_up_t);
		mdsd_dn.push(mdsd_dn_t);
		mdsd_bo.push(mdsd_bo_t);
	}
	conf.data.columns = [mdsd_up, mdsd_dn, mdsd_bo, mdsd_up_max];
	c3.generate(conf);
}
