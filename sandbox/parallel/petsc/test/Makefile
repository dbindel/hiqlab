PETSC_DIR=/usr/local/petsc/2.3.2-p6
PETSC_ARCH=linux_mpi1_cxx_hypre

NP               = 1

PROM_OPTS=\
	-pc_type prometheus \
	-pc_mg_type multiplicative -pc_mg_cycles 1 \
	-prometheus_preduce_base 500 \
	-prometheus_top_grid_limit 2500 \
	-aggmg_avoid_resmooth \
	-aggmg_smooths 1 \
	-prometheus_repartition \
	-malloc_dump

#	-info -mat_view_info
#	-on_error_attach_debugger
#	-info

include ${PETSC_DIR}/bmake/common/base

#-ksp_smonitor 

runhiq:
	-@${MPIRUN} -np 1 ../hiqlab -f dsb.lua \
	-ksp_gmres_cgs_refinement_type refine_always \
	${PROM_OPTS}

runhiqv:
	-@${MPIRUN} -np 1 valgrind --db-attach=yes --leak-check=full \
	../hiqlab -f dsb.lua -ksp_smonitor \
	-ksp_gmres_cgs_refinement_type refine_always ${PROM_OPTS}
