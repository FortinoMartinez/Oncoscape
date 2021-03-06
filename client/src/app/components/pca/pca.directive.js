(function() {
    'use strict';

    angular
        .module('oncoscape')
        .directive('osPca', explore);

    /** @ngInject */
    function explore() {

        var directive = {
            restrict: 'E',
            templateUrl: 'app/components/pca/pca.html',
            controller: PcaController,
            controllerAs: 'vm',
            bindToController: true
        };

        return directive;

        /** @ngInject */
        function PcaController($q, osApi, osHistory, $state, $stateParams, $timeout, $scope, d3, moment, $window, _) {

            // History
            var selectedIds = (osHistory.getPatientSelection() == null) ? [] : osHistory.getPatientSelection().ids;

            function saveSelected() {
                var selected = d3Chart.selectAll(".pca-node-selected")[0];
                if (selected.length == 0) return;
                osHistory.addPatientSelection("PCA", "Manual Selection",
                    d3Chart.selectAll(".pca-node-selected")[0].map(function(node) {
                        return node.__data__.id.toUpperCase();
                    })
                );
            }

            function setSelected() {
                if (selectedIds.length == 0) {
                    d3Chart.selectAll(".pca-node-selected").classed("pca-node-selected", false);
                } else {
                    d3Chart.selectAll("circle").classed("pca-node-selected", function() {
                        return (selectedIds.indexOf(this.__data__.id) >= 0)
                    });
                }
            }

            // Elements
            var d3Chart = d3.select("#pca-chart").append("svg").attr("id", "chart");
            var d3xAxis = d3Chart.append("g");
            var d3yAxis = d3Chart.append("g");
            var d3Tooltip = d3.select("body").append("div").attr("class", "tooltip pca-tooltip");

            // Properties
            var data;
            var layout = {
                width: 0,
                height: 0,
                xScale: 0,
                yScale: 0,
                xMax: 0,
                yMax: 0,
                xAxis: 0,
                yAxis: 0
            };

            // View Model
            var vm = (function(vm, osApi) {
                vm.datasource = osApi.getDataSource();
                vm.geneSets = [];
                vm.geneSet = null;
                vm.search = "";
                osApi.query("render_pca", {
                        disease: vm.datasource.disease,
                        $fields: ['geneset']
                    })
                    .then(function(response) {
                        vm.geneSets = response.data;
                        vm.geneSet = vm.geneSets[0];
                    });
                return vm;

            })(this, osApi)

            $scope.$watch('vm.geneSet', function(geneset) {
                if (geneset == null) return;
                osApi.query("render_pca", {
                        disease: vm.datasource.disease,
                        geneset: geneset.geneset
                    })
                    .then(function(response) {
                        vm.pc1 = response.data[0].pc1;
                        vm.pc2 = response.data[0].pc2;
                        var keys = Object.keys(response.data[0].data);
                        data = keys.map(function(key) {
                            this.data[key].push(key);
                            return this.data[key];
                        }, {
                            data: response.data[0].data
                        });
                        draw();
                    });
            });

            // Drawing Functions
            function scale() {
                layout.width = $window.innerWidth - 400;
                if (angular.element(".tray-right").attr("locked") == "false") {
                    layout.width += 300;
                }
                layout.height = $window.innerHeight - 190;
                if (angular.element(".tray").attr("locked") == "true") layout.width -= 300;

                d3Chart
                    .attr("width", '100%')
                    .attr("height", layout.height);
                layout.xScale = d3.scale.linear()
                    .domain([-layout.xMax, layout.xMax])
                    .range([0, layout.width]).nice();

                layout.yScale = d3.scale.linear()
                    .domain([-layout.yMax, layout.yMax])
                    .range([layout.height, 0]).nice();
            }

            function draw() {
                var vals = Object.keys(data).map(function(key) {
                    return data[key]
                }, {
                    data: data
                });
                layout.max = Math.abs(d3.max(vals, function(d) {
                    return +d[0];
                }));
                layout.min = Math.abs(d3.min(vals, function(d) {
                    return +d[0];
                }));
                layout.xMax = ((layout.max > layout.min) ? layout.max : layout.min) * 1.2;
                layout.max = Math.abs(d3.max(vals, function(d) {
                    return +d[1];
                }));
                layout.min = Math.abs(d3.min(vals, function(d) {
                    return +d[1];
                }));
                layout.yMax = ((layout.max > layout.min) ? layout.max : layout.min) * 1.2;

                // Refresh Scale
                scale();

                layout.xAxis = d3.svg.axis()
                    .scale(layout.xScale)
                    .orient("top")
                    .ticks(5);

                layout.yAxis = d3.svg.axis()
                    .scale(layout.yScale)
                    .orient("left")
                    .ticks(5);

                // Brush
                var brush = d3.svg.brush()
                    .x(layout.xScale)
                    .y(layout.yScale)
                    .on("brushend", function() {
                        var bv = brush.extent();
                        d3Chart.selectAll("circle")
                            .classed("pca-node-selected", function(d) {
                                return (d[0] > bv[0][0] && d[0] < bv[1][0] && d[1] > bv[0][1] && d[1] < bv[1][1]);
                            });
                        d3.select(this).transition().duration(300)
                            .call(brush.extent([
                                [0, 0],
                                [0, 0]
                            ]));
                        saveSelected();
                    });

                d3Chart.call(brush);

                var circles = d3Chart.selectAll("circle").data(data, function(d) {
                    return d;
                });

                circles.enter()
                    .append("circle")
                    .attr({
                        "class": "pca-node",
                        "cx": layout.width * .5,
                        "cy": layout.height * .5,
                        "r": 3
                    })
                    .style("fill-opacity", "0")
                    .on("mouseover", function(d) {
                        d3Tooltip.transition()
                            .duration(200)
                            .style("opacity", 1);
                        d3Tooltip.html(d.id)
                            .style("left", (d3.event.pageX + 10) + "px")
                            .style("top", (d3.event.pageY - 5) + "px");
                    })
                    .on("mouseout", function() {
                        d3Tooltip.transition()
                            .duration(500)
                            .style("opacity", 0);
                    })
                    .transition()
                    .duration(750)
                    .delay(function(d, i) {
                        return i / 300 * 500;
                    })
                    .attr("cx", function(d) {
                        return layout.xScale(d[0]);
                    })
                    .attr("cy", function(d) {
                        return layout.yScale(d[1]);
                    })
                    .style("fill-opacity", 1);

                circles.exit()
                    .transition()
                    .duration(600)
                    .delay(function(d, i) {
                        return i / 300 * 500;
                    })
                    .attr("cx", layout.width * .5)
                    .attr("cy", layout.height * .5)
                    .style("fill-opacity", "0")
                    .remove();


                d3yAxis
                    .attr("class", "axis")
                    .attr("transform", "translate(0, " + layout.yScale(0) + ")")
                    .call(layout.xAxis)
                    .append("text")
                    .text("PC1");

                d3xAxis
                    .attr("class", "axis")
                    .attr("transform", "translate(" + layout.xScale(0) + ", 0)")
                    .call(layout.yAxis)
                    .append("text")
                    .attr("y", 10)
                    .attr("dy", ".71em")
                    .text("PC2");

                setSelected();


            }


            vm.resize = function() {
                setScale();
                xAxis.scale(layout.xScale);
                yAxis.scale(layout.yScale);
                d3yAxis.attr("transform", "translate(0, " + yScale(0) + ")").call(layout.xAxis);
                d3xAxis.attr("transform", "translate(" + xScale(0) + ", 0)").call(layout.yAxis);
                d3Chart.selectAll("circle")
                    .attr("cx", function(d) {
                        return layout.xScale(d[0]);
                    })
                    .attr("cy", function(d) {
                        return layout.yScale(d[1]);
                    })
            };

            // Listen For Resize
            angular.element($window).bind('resize',
                _.debounce(vm.resize, 300)
            );
        }
    }
})();
