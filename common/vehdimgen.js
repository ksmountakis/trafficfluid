/* a,b,c,-segment according to EU classification

	link: https://www.automobiledimension.com/small-cars.php

	a-segment: avg length = 3.20m,	avg width = 1.6m,		ratio = 2.0
	b-segment: avg length = 3.90m,	avg width = 1.7m,		ratio = 2.3
	c-segment: avg length = 4.25m,	avg width = 1.75m,	ratio = 2.4
	d-segment: avg length = 4.55m,	avg width = 1.82m,	ratio = 2.5
	mid suv  : avg length = 4.60m,   avg width = 1.85m,	ratio = 2.5
	pickup tr: avg length = 5.15m,	avg width = 1.85m,	ratio = 2.8 
*/


// Relies on JStat.js

const MAXCARWID = 2;

function vehdimgen(type) {
	if (type === undefined) type = "car";

	if (type === "car") {
		var cases = [
			{len: 3.20, ratio: 2.0},
			{len: 3.90, ratio: 2.3},
			{len: 4.25, ratio: 2.4},
			{len: 4.55, ratio: 2.5},
			{len: 4.60, ratio: 2.5},
			{len: 5.15, ratio: 2.8}
		];
		let c = cases[Math.round(cases.length * jStat.beta.sample(3, 5))]; // mostly smaller cars
		let l = jStat.normal.sample(c.len, 0.1*c.len); // within [0.7*c.len, 1.3*c.len]
		let w = l / c.ratio; // always use average ratio

		return {l: l, w: w};
	}
}

