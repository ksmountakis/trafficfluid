let conf = 
{
	el: "#app",
	data: {
		sim: window.sim,
		t: 0,
		px_per_meter: 1,
		px_per_meter_max: 10,
		px_per_meter_min: 1,
		d3: {
			data: null,
            data_crashes: null,
            data_progress: 0
		},
        d3_data_progress: 0,
		play: {
			t: 0,
			active: false,
			intvl: null,
			speed: 2
		},
		 
		chart: {
			mode: "general",
			vehicle: {
				id: 0
			},
			general: {
				type: "mdsd"
			}
		}
	},

	methods: {
        road_container_resize(el) {
            let width_px = $(el).width() - 10;
            this.px_per_meter = Math.min(this.px_per_meter_max, Math.max(this.px_per_meter_min, width_px / this.sim.roadlen_meters));
            if (!this.d3.svg_main)
                this.svg_initialize();
            this.svg_show(this.play.t);
        },
        svg_initialize() {
            this.d3.svg_main = d3.select("#svg-main");
            this.d3.g_main = this.d3.svg_main.append("g")
                .attr("id", "g-main")
                .attr("transform", 
                `translate(0, 50) scale(${this.px_per_meter}, ${this.px_per_meter})`);

            // create road rectangle
            this.d3.road_rect = this.d3.g_main.append("rect")
                .attr("stroke", "white")
                .attr("stroke-width", 1.5/this.px_per_meter)
                .attr("x", "0").attr("y", "0")
                .attr("width", this.sim.roadlen_meters)
                .attr("height", this.sim.roadwid_meters)
                .attr("fill", "deepksyblue");
            
            // create frame boundaries
            this.d3.wall_y_up = this.d3.g_main.append("line")
                .attr("x1", 0).attr("x2", this.sim.roadlen_meters)
                .attr("y1", 0).attr("y2", 0)
                .attr("stroke", "skyblue")
                .attr("stroke-width", 1/this.px_per_meter);
            
            this.d3.wall_y_dn = this.d3.g_main.append("line")
                .attr("x1", 0).attr("x2", this.sim.roadlen_meters)
                .attr("y1", this.sim.roadwid_meters).attr("y2", this.sim.roadwid_meters)
                .attr("stroke", "skyblue")
                .attr("stroke-width", 1/this.px_per_meter);
            
            this.d3.wall_x = this.d3.g_main.append("line")
                .attr("x1", this.sim.roadlen_meters).attr("x2", this.sim.roadlen_meters)
                .attr("y1", 0).attr("y2", this.sim.roadwid_meters)
                .attr("stroke", "skyblue")
                .attr("stroke-width", 1/this.px_per_meter);

            this.svg_load_data();
            this.svg_show(0);
        },
        svg_load_data() {
            let pathmaker = d3.line().curve(d3.curveLinear);

            // create d3 data vector
            this.d3.data = [];
            this.d3.data_crashes = [];
            for (var t=0; t < this.sim.K; t++) {
                this.d3.data_progress = Math.round(((t/this.sim.K) * 100));
                this.d3.data.push([]);
                this.d3.data_crashes.push((t > 0)?this.d3.data_crashes[t-1]: 0);
                var j = 0;
                for (var i=0; i < this.sim.n; i++) {
                    // create data[t][i] = {}
                    this.d3.data[t].push({}); 
                    this.d3.data[t][j].id = j;
                    this.d3.data[t][j].l = this.sim.l[j];
                    this.d3.data[t][j].w = this.sim.w[j];
                    this.d3.data[t][j].vd = this.sim.vd[j];
                    this.d3.data[t][j].x = this.sim.x[j][t];
                    this.d3.data[t][j].y = this.sim.y[j][t];
                    this.d3.data[t][j].vx = this.sim.vx[j][t];
                    this.d3.data[t][j].crash = this.sim.crash[j][t];
                    if (this.sim.degree)
                        this.d3.data[t][j].degree = this.sim.degree[j][t];
                    this.d3.data_crashes[t] += this.d3.data[t][j].crash;

                    let leader = this.sim.leader[j][t];
                    var leader_points = [[0,0]];
                    this.d3.data[t][j].leader = leader;
                    if (leader !== -1 && this.sim.x[leader][t] > this.sim.x[j][t]) {
                        let arlen = 0.5;
                        let dx = this.sim.x[leader][t] -this.sim.x[j][t] - 1.1*arlen*2;
                        let dy = this.sim.y[leader][t] -this.sim.y[j][t] + (this.sim.w[leader])*0.5;
                    
                        leader_points = [
                            [this.sim.l[j], this.sim.w[j]*0.5], 
                            [dx, dy],
                            [dx, dy-arlen],
                            [dx+arlen*2, dy],
                            [dx, dy+arlen],
                            [dx, dy]
                        ];
                    }
                    
                    this.d3.data[t][j].leader_path = pathmaker(leader_points);
                    

                    this.d3.data[t][j].vx_mean = (t > 0)?this.d3.data[t-1][j].vx_mean + this.d3.data[t][j].vx: this.d3.data[0][j].vx;
                    if (this.sim.x[j][t] <= 2*this.play.speed*this.sim.vd[j]*this.sim.T || this.sim.x[j][t] >= this.sim.roadlen_meters-2*this.play.speed*this.sim.vd[j]*this.sim.T)
                        this.d3.data[t][j].reset = true;
                    else
                        this.d3.data[t][j].reset = false;
                    j++;
                }
            }

            for (var t=0; t < this.sim.K; t++) {
                for (var i=0; i < this.sim.n; i++) {
                    this.d3.data[t][i].vx_mean /= (t+1);
                }
            }
            this.d3.data_progress = 100;
        },
        svg_show(t) {
            const update = this.d3.g_main
                .selectAll("g.g-veh")
                .data(this.d3.data[t]);

            /* 
             * enter is non-empty only at svg_show(0), in which case:
             */

            /* a) we create a <g> associated with each datum in this.d3.data[]  */
            const enter = update.enter()
                    .append("g")
                    .attr("class", "g-veh");

            /* b) we create a <rect> within each such <g> */
            enter.append("rect")
                    .attr("stroke", "hotpink")
                    .attr("rx", 1/this.px_per_meter)
                    .attr("ry", 1/this.px_per_meter)
                    .attr("x", 0).attr("y", 0)
                    .attr("width", d => d.l)
                    .attr("height", d => d.w)
			  			  .on("click", (d, j) => {
							  if (this.chart.vehicle.id === j)
								  this.chart.vehicle.id = null;
							  else 
                                this.chart.vehicle.id = j;
                          });

            /* c) we create a <path> within each such <g> */
            /* */
            enter.append("path")
                .attr("fill", "grey")
                .attr("stroke", "grey")
                .attr("stroke-width", 1.5/this.px_per_meter);
            
            if (this.sim.degree) {
                enter.append("text")
                    .attr("font-size", 8/this.px_per_meter)
                    .attr("y", 7/this.px_per_meter)
                    .attr("x", 0/this.px_per_meter)
                    .attr("font-weight", "bold")
                    .attr("font-family", "monospace")
                    .attr("fill", "white")
                    .attr("class", "g-veh-text")
                    .text((d) => d.id);
            }
            /*
             * update.merge(enter) is non-empty for every svg_show(t) 
             */

            /* for each <g>, set vehicle position at t; note: (d.x, d.y) is center position */ 
            update.merge(enter)
                .attr("transform", d => `translate(${d.x-0.5*d.l}, ${d.y-0.5*d.w})`);

            /* for each <g> <rect>, set color according to deviation from desired speed */ 
            update.merge(enter)
                .select("rect")
                .attr("fill", (d,i) => {
                    if (i === this.chart.vehicle.id)
					    return "green";

                    if (i == 0 && sim.lanedrop !== 0 || d.crash !== 0)
                        return "black";
                    if (d.vx < 0.95*d.vd) 
                        return chroma.temperature(Math.pow(d.vx/(1+d.vd), 2) * 2500);
                    
                    if (d.vx <= 1.15*d.vd)
                        return "deepskyblue";

                    if (d.vx <= 1.35*d.vd)
                        return "blue";

                    return "purple";
                })
                .attr("stroke-width", (d) => (d.crash !== 0)?2/this.px_per_meter: 0);
    
            /* for each <g> <path>, set color according to deviation from desired speed */
            update.merge(enter)
                .select("path")
                .attr("d", (d) => d.leader_path);

            let j = this.chart.vehicle.id;
            if (j !== null) {
                
                this.d3.wall_x
                    .attr("x1", this.sim.wall_x[j][t])
                    .attr("x2", this.sim.wall_x[j][t])
                    .attr("y1", this.sim.wall_y_dn[j][t])
                    .attr("y2", this.sim.wall_y_up[j][t]);
                
                this.d3.wall_y_dn
                    .attr("x1", this.sim.x[j][t]-0.5*this.sim.l[j])
                    .attr("x2", this.sim.wall_x[j][t])
                    .attr("y1", this.sim.wall_y_dn[j][t])
                    .attr("y2", this.sim.wall_y_dn[j][t]);

                
                this.d3.wall_y_up
                    .attr("x1", this.sim.x[j][t]-0.5*this.sim.l[j])
                    .attr("x2", this.sim.wall_x[j][t])
                    .attr("y1", this.sim.wall_y_up[j][t])
                    .attr("y2", this.sim.wall_y_up[j][t]);
                
                if (this.sim.degree) {
                    update.merge(enter)
                        .select("text")
                        .text((d) => `${this.sim.degree[j][t][d.id].x},${this.sim.degree[j][t][d.id].y}`);
                }
            }

            update.exit().remove();
        },
        play_start() {
            this.play.active = true;
            this.svg_show(this.play.t);
            this.play.intvl = setInterval(() => {
                if (this.play.t >= this.sim.K-1) {
                    this.play_stop();
                    this.play.t = 0;
                    return;
                }
                this.play.t += 1;
                this.svg_show(this.play.t);
            }, this.sim.T*1000/this.play.speed);
        },
        play_stop() {
            this.play.active = false;
            clearInterval(this.play.intvl);
        },
        chart_show() {
            window.chart_show(this);
        }
	}
};

Vue.use(Vuetify);
new Vue(conf);
