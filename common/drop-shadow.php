<?php // vim: set filetype=html : ?>

<!-- drop shadow filter -->
<svg>
	<filter id="drop-shadow"  width="150%" height="150%">
		<feGaussianBlur in="SourceAlpha" stdDeviation="2.2"/>
		<feOffset dx="2" dy="2" result="offsetblur"/>
		<feFlood flood-color="rgba(0,0,0,0.5)"/>
		<feComposite in2="offsetblur" operator="in"/>
		<feMerge>
			<feMergeNode/>
			<feMergeNode in="SourceGraphic"/>
		</feMerge>
	</filter>
</svg>

<style>
	.drop-shadow {
		filter: url(#drop-shadow);
	}
	.fade-edges {
		box-shadow: 0 0 8px 8px white inset;
	}
</style>
